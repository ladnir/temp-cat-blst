




#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/fields/Fp12.h"
#include <string.h>
#include "blst/fields/Fr.h"

#ifndef BLST_FP12_TOWER_H
#define BLST_FP12_TOWER_H

#ifndef FP12_CPP_INLINE
#define FP12_CPP_INLINE
#endif // !FP12_CPP_INLINE


namespace blst
{




#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
#define __FP2x2__
	/*
	 * Fp2x2 is a "widened" version of Fp2, which allows to consolidate
	 * reductions from several multiplications. In other words instead of
	 * "mul_redc-mul_redc-add" we get "mul-mul-add-redc," where latter
	 * addition is double-width... To be more specific this gives ~7-10%
	 * faster pairing depending on platform...
	 */

	BLST_HD BLST_INLINE void add_fp2x2(UFp2& ret, const UFp2& a, const UFp2& b)
	{
		add_mod_384x384(&ret[0], &a[0], &b[0], &get_BLST_P());
		add_mod_384x384(&ret[1], &a[1], &b[1], &get_BLST_P());
	}

	BLST_HD BLST_INLINE void sub_fp2x2(UFp2& ret, const UFp2& a, const UFp2& b)
	{
		sub_mod_384x384(&ret[0], &a[0], &b[0], &get_BLST_P());
		sub_mod_384x384(&ret[1], &a[1], &b[1], &get_BLST_P());
	}

	BLST_HD BLST_INLINE void mul_by_u_plus_1_fp2x2(UFp2& ret, const UFp2& a)
	{
		/* caveat lector! |ret| may not be same as |a| */
		sub_mod_384x384(&ret[0], &a[0], &a[1], &get_BLST_P());
		add_mod_384x384(&ret[1], &a[0], &a[1], &get_BLST_P());
	}

	BLST_HD BLST_INLINE void redc_fp2x2(Fp2& ret, const UFp2& a)
	{
		Fp::redc(ret[0], a[0]);
		Fp::redc(ret[1], a[1]);
	}

	BLST_HD BLST_INLINE void mul_fp2x2(UFp2& ret, const Fp2& a, const Fp2& b)
	{
#if 1
		mul_v382x(ret, a, b, get_BLST_P());   /* +~6% in Miller loop */
#else
		union { Fp x[2]; UFp& x2; } t;

		add_mod_384(t.x[0], a[0], a[1], get_BLST_P());
		add_mod_384(t.x[1], b[0], b[1], get_BLST_P());
		mul_384(ret[1], t.x[0], t.x[1]);

		mul_384(ret[0], a[0], b[0]);
		mul_384(t.x2, a[1], b[1]);

		sub_mod_384x384(ret[1], ret[1], ret[0], get_BLST_P());
		sub_mod_384x384(ret[1], ret[1], t.x2, get_BLST_P());

		sub_mod_384x384(ret[0], ret[0], t.x2, get_BLST_P());
#endif
	}

	BLST_HD BLST_INLINE void sqr_fp2x2(UFp2& ret, const Fp2& a)
	{
#if 1
		sqr_v382x(ret, a, get_BLST_P());      /* +~5% in final exponentiation */
#else
		Fp& t0, t1;

		add_mod_384(t0, a[0], a[1], get_BLST_P());
		sub_mod_384(t1, a[0], a[1], get_BLST_P());

		mul_384(ret[1], a[0], a[1]);
		add_mod_384x384(ret[1], ret[1], ret[1], get_BLST_P());

		mul_384(ret[0], t0, t1);
#endif
	}

	//typedef UFp2 UFp6&[3];
	struct UFp6 {
		UFp2 data[3];

		UFp2& operator[](size_t idx) {
			return data[idx];
		}
		const UFp2& operator[](size_t idx) const {
			return data[idx];
		}
	};

	BLST_HD BLST_INLINE void sub_fp6x2(UFp6& ret, const UFp6& a,
		const UFp6& b)
	{
		sub_fp2x2(ret[0], a[0], b[0]);
		sub_fp2x2(ret[1], a[1], b[1]);
		sub_fp2x2(ret[2], a[2], b[2]);
	}

	BLST_HD BLST_INLINE void mul_fp6x2(UFp6& ret, const Fp6& a, const Fp6& b)
	{
		UFp2 t0, t1, t2;
		Fp2 aa, bb;
		const Fp2& a0 = a[0];
		mul_fp2x2(t0, a0, b[0]);
		mul_fp2x2(t1, a[1], b[1]);
		mul_fp2x2(t2, a[2], b[2]);

		/* ret[0] = ((a1 + a2)*(b1 + b2) - a1*b1 - a2*b2)*(u+1) + a0*b0
				  = (a1*b2 + a2*b1)*(u+1) + a0*b0 */
		Fp2::add(aa, a[1], a[2]);
		Fp2::add(bb, b[1], b[2]);
		mul_fp2x2(ret[0], aa, bb);
		sub_fp2x2(ret[0], ret[0], t1);
		sub_fp2x2(ret[0], ret[0], t2);
		mul_by_u_plus_1_fp2x2(ret[1], ret[0]);  /* borrow ret[1] for a moment */
		add_fp2x2(ret[0], ret[1], t0);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*b2*(u+1)
				  = a0*b1 + a1*b0 + a2*b2*(u+1) */
		Fp2::add(aa, a[0], a[1]);
		Fp2::add(bb, b[0], b[1]);
		mul_fp2x2(ret[1], aa, bb);
		sub_fp2x2(ret[1], ret[1], t0);
		sub_fp2x2(ret[1], ret[1], t1);
		mul_by_u_plus_1_fp2x2(ret[2], t2);      /* borrow ret[2] for a moment */
		add_fp2x2(ret[1], ret[1], ret[2]);

		/* ret[2] = (a0 + a2)*(b0 + b2) - a0*b0 - a2*b2 + a1*b1
				  = a0*b2 + a2*b0 + a1*b1 */
		Fp2::add(aa, a[0], a[2]);
		Fp2::add(bb, b[0], b[2]);
		mul_fp2x2(ret[2], aa, bb);
		sub_fp2x2(ret[2], ret[2], t0);
		sub_fp2x2(ret[2], ret[2], t2);
		add_fp2x2(ret[2], ret[2], t1);
	}

	BLST_HD BLST_INLINE void redc_fp6x2(Fp6& ret, const UFp6& a)
	{
		redc_fp2x2(ret[0], a[0]);
		redc_fp2x2(ret[1], a[1]);
		redc_fp2x2(ret[2], a[2]);
	}

	BLST_HD BLST_INLINE void mul_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
	{
		UFp6 r;

		mul_fp6x2(r, a, b);
		redc_fp6x2(ret, r); /* narrow to normal width */
	}

	BLST_HD BLST_INLINE void sqr_fp6(Fp6& ret, const Fp6& a)
	{
		UFp2 s0, m01, m12, s2, rx;

		sqr_fp2x2(s0, a[0]);

		mul_fp2x2(m01, a[0], a[1]);
		add_fp2x2(m01, m01, m01);

		mul_fp2x2(m12, a[1], a[2]);
		add_fp2x2(m12, m12, m12);

		sqr_fp2x2(s2, a[2]);

		/* ret[2] = (a0 + a1 + a2)^2 - a0^2 - a2^2 - 2*(a0*a1) - 2*(a1*a2)
				  = a1^2 + 2*(a0*a2) */
		Fp2::add(ret[2], a[2], a[1]);
		Fp2::add(ret[2], ret[2], a[0]);
		sqr_fp2x2(rx, ret[2]);
		sub_fp2x2(rx, rx, s0);
		sub_fp2x2(rx, rx, s2);
		sub_fp2x2(rx, rx, m01);
		sub_fp2x2(rx, rx, m12);
		redc_fp2x2(ret[2], rx);

		/* ret[0] = a0^2 + 2*(a1*a2)*(u+1) */
		mul_by_u_plus_1_fp2x2(rx, m12);
		add_fp2x2(rx, rx, s0);
		redc_fp2x2(ret[0], rx);

		/* ret[1] = a2^2*(u+1) + 2*(a0*a1) */
		mul_by_u_plus_1_fp2x2(rx, s2);
		add_fp2x2(rx, rx, m01);
		redc_fp2x2(ret[1], rx);
	}



#else // __FP2x2__






	BLST_HD BLST_INLINE void mul_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
	{
		Fp2 t0, t1, t2, t3, t4, t5;

		Fp2::mul(t0, a[0], b[0]);
		Fp2::mul(t1, a[1], b[1]);
		Fp2::mul(t2, a[2], b[2]);

		/* ret[0] = ((a1 + a2)*(b1 + b2) - a1*b1 - a2*b2)*(u+1) + a0*b0
				  = (a1*b2 + a2*b1)*(u+1) + a0*b0 */
		Fp2::add(t4, a[1], a[2]);
		Fp2::add(t5, b[1], b[2]);
		Fp2::mul(t3, t4, t5);
		Fp2::sub(t3, t3, t1);
		Fp2::sub(t3, t3, t2);
		Fp2::mul_by_u_plus_1(t3, t3);
		/* Fp2::add(ret[0], t3, t0); considering possible aliasing... */

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*b2*(u+1)
				  = a0*b1 + a1*b0 + a2*b2*(u+1) */
		Fp2::add(t4, a[0], a[1]);
		Fp2::add(t5, b[0], b[1]);
		Fp2::mul(ret[1], t4, t5);
		Fp2::sub(ret[1], ret[1], t0);
		Fp2::sub(ret[1], ret[1], t1);
		Fp2::mul_by_u_plus_1(t4, t2);
		Fp2::add(ret[1], ret[1], t4);

		/* ret[2] = (a0 + a2)*(b0 + b2) - a0*b0 - a2*b2 + a1*b1
				  = a0*b2 + a2*b0 + a1*b1 */
		Fp2::add(t4, a[0], a[2]);
		Fp2::add(t5, b[0], b[2]);
		Fp2::mul(ret[2], t4, t5);
		Fp2::sub(ret[2], ret[2], t0);
		Fp2::sub(ret[2], ret[2], t2);
		Fp2::add(ret[2], ret[2], t1);

		Fp2::add(ret[0], t3, t0);    /* ... moved from above */
	}

	BLST_HD BLST_INLINE void sqr_fp6(Fp6& ret, const Fp6& a)
	{
		Fp2 s0, m01, m12, s2;

		Fp2::sqr(s0, a[0]);

		Fp2::mul(m01, a[0], a[1]);
		Fp2::add(m01, m01, m01);

		Fp2::mul(m12, a[1], a[2]);
		Fp2::add(m12, m12, m12);

		Fp2::sqr(s2, a[2]);

		/* ret[2] = (a0 + a1 + a2)^2 - a0^2 - a2^2 - 2*(a0*a1) - 2*(a1*a2)
				  = a1^2 + 2*(a0*a2) */
		Fp2::add(ret[2], a[2], a[1]);
		Fp2::add(ret[2], ret[2], a[0]);
		Fp2::sqr(ret[2], ret[2]);
		Fp2::sub(ret[2], ret[2], s0);
		Fp2::sub(ret[2], ret[2], s2);
		Fp2::sub(ret[2], ret[2], m01);
		Fp2::sub(ret[2], ret[2], m12);

		/* ret[0] = a0^2 + 2*(a1*a2)*(u+1) */
		Fp2::mul_by_u_plus_1(ret[0], m12);
		Fp2::add(ret[0], ret[0], s0);

		/* ret[1] = a2^2*(u+1) + 2*(a0*a1) */
		Fp2::mul_by_u_plus_1(ret[1], s2);
		Fp2::add(ret[1], ret[1], m01);
	}
#endif // __FP2x2__

	BLST_HD BLST_INLINE void add_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
	{
		Fp2::add(ret[0], a[0], b[0]);
		Fp2::add(ret[1], a[1], b[1]);
		Fp2::add(ret[2], a[2], b[2]);
	}

	BLST_HD BLST_INLINE void sub_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
	{
		Fp2::sub(ret[0], a[0], b[0]);
		Fp2::sub(ret[1], a[1], b[1]);
		Fp2::sub(ret[2], a[2], b[2]);
	}

	BLST_HD BLST_INLINE void neg_fp6(Fp6& ret, const Fp6& a)
	{
		Fp2::neg(ret[0], a[0]);
		Fp2::neg(ret[1], a[1]);
		Fp2::neg(ret[2], a[2]);
	}




	BLST_HD BLST_INLINE void Fp12::add(Fp12& ret, const Fp12& a, const Fp12& b)
	{
		add_fp6(ret[0], a[0], b[0]);
		add_fp6(ret[1], a[1], b[1]);
	}
	BLST_HD BLST_INLINE void Fp12::sub(Fp12& ret, const Fp12& a, const Fp12& b)
	{
		sub_fp6(ret[0], a[0], b[0]);
		sub_fp6(ret[1], a[1], b[1]);
	}


	//#if 0
	//#define mul_by_v_fp6 mul_by_v_fp6
	//    BLST_HD BLST_INLINE void mul_by_v_fp6(Fp6& ret, const Fp6& a)
	//    {
	//        Fp2 t;
	//
	//        Fp2::mul_by_u_plus_1(t, a[2]);
	//        vec_copy(ret[2], a[1], sizeof(a[1]));
	//        vec_copy(ret[1], a[0], sizeof(a[0]));
	//        vec_copy(ret[0], t, sizeof(t));
	//    }
	//#endif

		/*
		 * Fp12 extension
		 */
#if defined(__FP2x2__)
	BLST_HD FP12_CPP_INLINE void Fp12::mul(Fp12& ret, const Fp12& a, const Fp12& b)
	{
		UFp6 t0, t1, rx;
		Fp6 t2;

		mul_fp6x2(t0, a[0], b[0]);
		mul_fp6x2(t1, a[1], b[1]);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
				  = a0*b1 + a1*b0 */
		add_fp6(t2, a[0], a[1]);
		add_fp6(ret[1], b[0], b[1]);
		mul_fp6x2(rx, ret[1], t2);
		sub_fp6x2(rx, rx, t0);
		sub_fp6x2(rx, rx, t1);
		redc_fp6x2(ret[1], rx);

		/* ret[0] = a0*b0 + a1*b1*v */
		mul_by_u_plus_1_fp2x2(rx[0], t1[2]);
		add_fp2x2(rx[0], t0[0], rx[0]);
		add_fp2x2(rx[1], t0[1], t1[0]);
		add_fp2x2(rx[2], t0[2], t1[1]);
		redc_fp6x2(ret[0], rx);
	}

	BLST_HD BLST_INLINE void mul_by_0y0_fp6x2(UFp6& ret, const Fp6& a,
		const Fp2& b)
	{
		mul_fp2x2(ret[1], a[2], b);     /* borrow ret[1] for a moment */
		mul_by_u_plus_1_fp2x2(ret[0], ret[1]);
		mul_fp2x2(ret[1], a[0], b);
		mul_fp2x2(ret[2], a[1], b);
	}

	BLST_HD BLST_INLINE void mul_by_xy0_fp6x2(UFp6& ret, const Fp6& a,
		const Fp6& b)
	{
		UFp2 t0, t1;
		Fp2 aa, bb;

		mul_fp2x2(t0, a[0], b[0]);
		mul_fp2x2(t1, a[1], b[1]);

		/* ret[0] = ((a1 + a2)*(b1 + 0) - a1*b1 - a2*0)*(u+1) + a0*b0
				  = (a1*0 + a2*b1)*(u+1) + a0*b0 */
		mul_fp2x2(ret[1], a[2], b[1]);  /* borrow ret[1] for a moment */
		mul_by_u_plus_1_fp2x2(ret[0], ret[1]);
		add_fp2x2(ret[0], ret[0], t0);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*0*(u+1)
				  = a0*b1 + a1*b0 + a2*0*(u+1) */
		Fp2::add(aa, a[0], a[1]);
		Fp2::add(bb, b[0], b[1]);
		mul_fp2x2(ret[1], aa, bb);
		sub_fp2x2(ret[1], ret[1], t0);
		sub_fp2x2(ret[1], ret[1], t1);

		/* ret[2] = (a0 + a2)*(b0 + 0) - a0*b0 - a2*0 + a1*b1
				  = a0*0 + a2*b0 + a1*b1 */
		mul_fp2x2(ret[2], a[2], b[0]);
		add_fp2x2(ret[2], ret[2], t1);
	}

	BLST_HD FP12_CPP_INLINE void Fp12::mul_by_xy00z0(Fp12& ret, const Fp12& a,
		const Fp6& xy00z0)
	{
		UFp6 t0, t1, rr;
		Fp6 t2;

		mul_by_xy0_fp6x2(t0, a[0], xy00z0);
		mul_by_0y0_fp6x2(t1, a[1], xy00z0[2]);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
				  = a0*b1 + a1*b0 */
		vec_copy(&t2[0], &xy00z0[0], sizeof(t2[0]));
		Fp2::add(t2[1], xy00z0[1], xy00z0[2]);
		add_fp6(ret[1], a[0], a[1]);
		mul_by_xy0_fp6x2(rr, ret[1], t2);
		sub_fp6x2(rr, rr, t0);
		sub_fp6x2(rr, rr, t1);
		redc_fp6x2(ret[1], rr);

		/* ret[0] = a0*b0 + a1*b1*v */
		mul_by_u_plus_1_fp2x2(rr[0], t1[2]);
		add_fp2x2(rr[0], t0[0], rr[0]);
		add_fp2x2(rr[1], t0[1], t1[0]);
		add_fp2x2(rr[2], t0[2], t1[1]);
		redc_fp6x2(ret[0], rr);
	}
#else // __FP2x2__
	BLST_HD FP12_CPP_INLINE void Fp12::mul(Fp12& ret, const Fp12& a, const Fp12& b)
	{
		Fp6 t0, t1, t2;

		mul_fp6(t0, a[0], b[0]);
		mul_fp6(t1, a[1], b[1]);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
				  = a0*b1 + a1*b0 */
		add_fp6(t2, a[0], a[1]);
		add_fp6(ret[1], b[0], b[1]);
		mul_fp6(ret[1], ret[1], t2);
		sub_fp6(ret[1], ret[1], t0);
		sub_fp6(ret[1], ret[1], t1);

		/* ret[0] = a0*b0 + a1*b1*v */
#ifdef mul_by_v_fp6
		mul_by_v_fp6(t1, t1);
		add_fp6(ret[0], t0, t1);
#else
		Fp2::mul_by_u_plus_1(t1[2], t1[2]);
		Fp2::add(ret[0][0], t0[0], t1[2]);
		Fp2::add(ret[0][1], t0[1], t1[0]);
		Fp2::add(ret[0][2], t0[2], t1[1]);
#endif
	}

	BLST_HD FP12_CPP_INLINE void mul_by_0y0_fp6(Fp6& ret, const Fp6& a,
		const Fp2& b)
	{
		Fp2 t;

		Fp2::mul(t, a[2], b);
		Fp2::mul(ret[2], a[1], b);
		Fp2::mul(ret[1], a[0], b);
		Fp2::mul_by_u_plus_1(ret[0], t);
	}

	BLST_HD FP12_CPP_INLINE void mul_by_xy0_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
	{
		Fp2 t0, t1, /*t2,*/ t3, t4, t5;

		Fp2::mul(t0, a[0], b[0]);
		Fp2::mul(t1, a[1], b[1]);

		/* ret[0] = ((a1 + a2)*(b1 + 0) - a1*b1 - a2*0)*(u+1) + a0*b0
				  = (a1*0 + a2*b1)*(u+1) + a0*b0 */
		Fp2::mul(t3, a[2], b[1]);
		Fp2::mul_by_u_plus_1(t3, t3);
		/* Fp2::add(ret[0], t3, t0); considering possible aliasing... */

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*0*(u+1)
				  = a0*b1 + a1*b0 + a2*0*(u+1) */
		Fp2::add(t4, a[0], a[1]);
		Fp2::add(t5, b[0], b[1]);
		Fp2::mul(ret[1], t4, t5);
		Fp2::sub(ret[1], ret[1], t0);
		Fp2::sub(ret[1], ret[1], t1);

		/* ret[2] = (a0 + a2)*(b0 + 0) - a0*b0 - a2*0 + a1*b1
				  = a0*0 + a2*b0 + a1*b1 */
		Fp2::mul(ret[2], a[2], b[0]);
		Fp2::add(ret[2], ret[2], t1);

		Fp2::add(ret[0], t3, t0);    /* ... moved from above */
	}

	BLST_HD FP12_CPP_INLINE void Fp12::mul_by_xy00z0(Fp12& ret, const Fp12& a,
		const Fp6& xy00z0)
	{
		Fp6 t0, t1, t2;

		mul_by_xy0_fp6(t0, a[0], xy00z0);
		mul_by_0y0_fp6(t1, a[1], xy00z0[2]);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
				  = a0*b1 + a1*b0 */
		vec_copy(&t2[0], &xy00z0[0], sizeof(t2[0]));
		Fp2::add(t2[1], xy00z0[1], xy00z0[2]);
		add_fp6(ret[1], a[0], a[1]);
		mul_by_xy0_fp6(ret[1], ret[1], t2);
		sub_fp6(ret[1], ret[1], t0);
		sub_fp6(ret[1], ret[1], t1);

		/* ret[0] = a0*b0 + a1*b1*v */
#ifdef mul_by_v_fp6
		mul_by_v_fp6(t1, t1);
		add_fp6(ret[0], t0, t1);
#else
		Fp2::mul_by_u_plus_1(t1[2], t1[2]);
		Fp2::add(ret[0][0], t0[0], t1[2]);
		Fp2::add(ret[0][1], t0[1], t1[0]);
		Fp2::add(ret[0][2], t0[2], t1[1]);
#endif
	}
#endif // __FP2x2__


#if defined(__FP2x2__)
	
	// Helper for the FP2x2 optimization path
	BLST_HD BLST_INLINE void add_fp6x2(UFp6& ret, const UFp6& a, const UFp6& b)
	{
		add_fp2x2(ret[0], a[0], b[0]);
		add_fp2x2(ret[1], a[1], b[1]);
		add_fp2x2(ret[2], a[2], b[2]);
	}

	// Optimized mul with integrated conjugation
	// Handles both __FP2x2__ (CPU/ASM) and Standard (CUDA/Portable) paths
	BLST_HD FP12_CPP_INLINE void Fp12::mul(Fp12& ret, const Fp12& a, const Fp12& b, bool conj_b)
	{
		UFp6 t0, t1, rx;
		Fp6 t2;

		mul_fp6x2(t0, a[0], b[0]);
		mul_fp6x2(t1, a[1], b[1]);

		/* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 */
		add_fp6(t2, a[0], a[1]);

		// Optimization: Merge conjugation into the linear combination
		if (conj_b) sub_fp6(ret[1], b[0], b[1]); // b0 + (-b1)
		else        add_fp6(ret[1], b[0], b[1]);

		mul_fp6x2(rx, ret[1], t2);
		sub_fp6x2(rx, rx, t0);

		if (conj_b) add_fp6x2(rx, rx, t1); // - (a1 * -b1) = + a1*b1
		else        sub_fp6x2(rx, rx, t1);

		redc_fp6x2(ret[1], rx);

		/* ret[0] = a0*b0 +/- a1*b1*v */
		// Note: t1[2] corresponds to high part of expansion used for 'v' mul
		mul_by_u_plus_1_fp2x2(rx[0], t1[2]);

		if (conj_b) {
			// a1*b1*v term is subtracted
			sub_fp2x2(rx[0], t0[0], rx[0]);
			sub_fp2x2(rx[1], t0[1], t1[0]);
			sub_fp2x2(rx[2], t0[2], t1[1]);
		}
		else {
			add_fp2x2(rx[0], t0[0], rx[0]);
			add_fp2x2(rx[1], t0[1], t1[0]);
			add_fp2x2(rx[2], t0[2], t1[1]);
		}
		redc_fp6x2(ret[0], rx);
	}
#else // Normal path (Likely active for CUDA)
	BLST_HD FP12_CPP_INLINE void Fp12::mul(Fp12& ret, const Fp12& a, const Fp12& b, bool conj_b)
	{
		Fp6 t0, t1, t2;

		mul_fp6(t0, a[0], b[0]);
		mul_fp6(t1, a[1], b[1]);

		/* ret[1] calculation */
		add_fp6(t2, a[0], a[1]);

		if (conj_b) sub_fp6(ret[1], b[0], b[1]); // b0 - b1
		else        add_fp6(ret[1], b[0], b[1]);

		mul_fp6(ret[1], ret[1], t2);
		sub_fp6(ret[1], ret[1], t0);

		if (conj_b) add_fp6(ret[1], ret[1], t1); // - (a1 * -b1)
		else        sub_fp6(ret[1], ret[1], t1);

		/* ret[0] calculation */
#ifdef mul_by_v_fp6
		mul_by_v_fp6(t1, t1);
		if (conj_b) sub_fp6(ret[0], t0, t1);
		else        add_fp6(ret[0], t0, t1);
#else
		Fp2::mul_by_u_plus_1(t1[2], t1[2]);
		if (conj_b) {
			Fp2::sub(ret[0][0], t0[0], t1[2]);
			Fp2::sub(ret[0][1], t0[1], t1[0]);
			Fp2::sub(ret[0][2], t0[2], t1[1]);
		}
		else {
			Fp2::add(ret[0][0], t0[0], t1[2]);
			Fp2::add(ret[0][1], t0[1], t1[0]);
			Fp2::add(ret[0][2], t0[2], t1[1]);
		}
#endif
	}
#endif

	BLST_HD FP12_CPP_INLINE void Fp12::sqr(Fp12& ret, const Fp12& a)
	{
		Fp6 t0, t1;

		add_fp6(t0, a[0], a[1]);
#ifdef mul_by_v_fp6
		mul_by_v_fp6(t1, a[1]);
		add_fp6(t1, a[0], t1);
#else
		Fp2::mul_by_u_plus_1(t1[2], a[1][2]);
		Fp2::add(t1[0], a[0][0], t1[2]);
		Fp2::add(t1[1], a[0][1], a[1][0]);
		Fp2::add(t1[2], a[0][2], a[1][1]);
#endif
		mul_fp6(t0, t0, t1);
		mul_fp6(t1, a[0], a[1]);

		/* ret[1] = 2*(a0*a1) */
		add_fp6(ret[1], t1, t1);

		/* ret[0] = (a0 + a1)*(a0 + a1*v) - a0*a1 - a0*a1*v
				  = a0^2 + a1^2*v */
		sub_fp6(ret[0], t0, t1);
#ifdef mul_by_v_fp6
		mul_by_v_fp6(t1, t1);
		sub_fp6(ret[0], ret[0], t1);
#else
		Fp2::mul_by_u_plus_1(t1[2], t1[2]);
		Fp2::sub(ret[0][0], ret[0][0], t1[2]);
		Fp2::sub(ret[0][1], ret[0][1], t1[0]);
		Fp2::sub(ret[0][2], ret[0][2], t1[1]);
#endif
	}

	BLST_HD FP12_CPP_INLINE void Fp12::conjugate(Fp12& a)
	{
		neg_fp6(a[1], a[1]);
	}

	BLST_HD FP12_CPP_INLINE void inverse_fp6(Fp6& ret, const Fp6& a)
	{
		Fp2 c0, c1, c2, t0, t1;

		/* c0 = a0^2 - (a1*a2)*(u+1) */
		Fp2::sqr(c0, a[0]);
		Fp2::mul(t0, a[1], a[2]);
		Fp2::mul_by_u_plus_1(t0, t0);
		Fp2::sub(c0, c0, t0);

		/* c1 = a2^2*(u+1) - (a0*a1) */
		Fp2::sqr(c1, a[2]);
		Fp2::mul_by_u_plus_1(c1, c1);
		Fp2::mul(t0, a[0], a[1]);
		Fp2::sub(c1, c1, t0);

		/* c2 = a1^2 - a0*a2 */
		Fp2::sqr(c2, a[1]);
		Fp2::mul(t0, a[0], a[2]);
		Fp2::sub(c2, c2, t0);

		/* (a2*c1 + a1*c2)*(u+1) + a0*c0 */
		Fp2::mul(t0, c1, a[2]);
		Fp2::mul(t1, c2, a[1]);
		Fp2::add(t0, t0, t1);
		Fp2::mul_by_u_plus_1(t0, t0);
		Fp2::mul(t1, c0, a[0]);
		Fp2::add(t0, t0, t1);

		Fp2::reciprocal(t1, t0);

		Fp2::mul(ret[0], c0, t1);
		Fp2::mul(ret[1], c1, t1);
		Fp2::mul(ret[2], c2, t1);
	}

	BLST_HD FP12_CPP_INLINE void Fp12::inverse(Fp12& ret, const Fp12& a)
	{
		Fp6 t0, t1;

		sqr_fp6(t0, a[0]);
		sqr_fp6(t1, a[1]);
#ifdef mul_by_v_fp6
		mul_by_v_fp6(t1, t1);
		sub_fp6(t0, t0, t1);
#else
		Fp2::mul_by_u_plus_1(t1[2], t1[2]);
		Fp2::sub(t0[0], t0[0], t1[2]);
		Fp2::sub(t0[1], t0[1], t1[0]);
		Fp2::sub(t0[2], t0[2], t1[1]);
#endif

		inverse_fp6(t1, t0);

		mul_fp6(ret[0], a[0], t1);
		mul_fp6(ret[1], a[1], t1);
		neg_fp6(ret[1], ret[1]);
	}

	//typedef Fp2 Fpfp4[2];

#if defined(__FP2x2__)
	BLST_HD FP12_CPP_INLINE void sqr_fp4(Fp4& ret, const Fp2& a0, const Fp2& a1)
	{
		UFp2 t0, t1, t2;

		sqr_fp2x2(t0, a0);
		sqr_fp2x2(t1, a1);
		Fp2::add(ret[1], a0, a1);

		mul_by_u_plus_1_fp2x2(t2, t1);
		add_fp2x2(t2, t2, t0);
		redc_fp2x2(ret[0], t2);

		sqr_fp2x2(t2, ret[1]);
		sub_fp2x2(t2, t2, t0);
		sub_fp2x2(t2, t2, t1);
		redc_fp2x2(ret[1], t2);
	}
#else //__FP2x2__
	BLST_HD FP12_CPP_INLINE void sqr_fp4(Fp4& ret, const Fp2& a0, const Fp2& a1)
	{
		Fp2 t0, t1;

		Fp2::sqr(t0, a0);
		Fp2::sqr(t1, a1);
		Fp2::add(ret[1], a0, a1);

		Fp2::mul_by_u_plus_1(ret[0], t1);
		Fp2::add(ret[0], ret[0], t0);

		Fp2::sqr(ret[1], ret[1]);
		Fp2::sub(ret[1], ret[1], t0);
		Fp2::sub(ret[1], ret[1], t1);
	}
#endif // __FP2x2__

	BLST_HD FP12_CPP_INLINE void Fp12::cyclotomic_sqr(Fp12& ret, const Fp12& a)
	{
		if (0)
		{

			Fp4 t0, t1, t2;

			sqr_fp4(t0, a[0][0], a[1][1]);
			sqr_fp4(t1, a[1][0], a[0][2]);
			sqr_fp4(t2, a[0][1], a[1][2]);

			Fp2::sub(ret[0][0], t0[0], a[0][0]);
			Fp2::add(ret[0][0], ret[0][0], ret[0][0]);
			Fp2::add(ret[0][0], ret[0][0], t0[0]);

			Fp2::sub(ret[0][1], t1[0], a[0][1]);
			Fp2::add(ret[0][1], ret[0][1], ret[0][1]);
			Fp2::add(ret[0][1], ret[0][1], t1[0]);

			Fp2::sub(ret[0][2], t2[0], a[0][2]);
			Fp2::add(ret[0][2], ret[0][2], ret[0][2]);
			Fp2::add(ret[0][2], ret[0][2], t2[0]);

			Fp2::mul_by_u_plus_1(t2[1], t2[1]);
			Fp2::add(ret[1][0], t2[1], a[1][0]);
			Fp2::add(ret[1][0], ret[1][0], ret[1][0]);
			Fp2::add(ret[1][0], ret[1][0], t2[1]);

			Fp2::add(ret[1][1], t0[1], a[1][1]);
			Fp2::add(ret[1][1], ret[1][1], ret[1][1]);
			Fp2::add(ret[1][1], ret[1][1], t0[1]);

			Fp2::add(ret[1][2], t1[1], a[1][2]);
			Fp2::add(ret[1][2], ret[1][2], ret[1][2]);
			Fp2::add(ret[1][2], ret[1][2], t1[1]);
		}

		if (1)
		{
			{

				Fp4 t0;
				sqr_fp4(t0, a[0][0], a[1][1]);
				Fp2::sub(ret[0][0], t0[0], a[0][0]);
				Fp2::add(ret[0][0], ret[0][0], ret[0][0]);
				Fp2::add(ret[0][0], ret[0][0], t0[0]);

				Fp2::add(ret[1][1], t0[1], a[1][1]);
				Fp2::add(ret[1][1], ret[1][1], ret[1][1]);
				Fp2::add(ret[1][1], ret[1][1], t0[1]);

			}

			// must compute this one early in case &a==&ret.
			Fp4 t2;
			sqr_fp4(t2, a[0][1], a[1][2]);

			{
				Fp4 t1;
				sqr_fp4(t1, a[1][0], a[0][2]);
				Fp2::sub(ret[0][1], t1[0], a[0][1]);
				Fp2::add(ret[0][1], ret[0][1], ret[0][1]);
				Fp2::add(ret[0][1], ret[0][1], t1[0]);

				Fp2::add(ret[1][2], t1[1], a[1][2]);
				Fp2::add(ret[1][2], ret[1][2], ret[1][2]);
				Fp2::add(ret[1][2], ret[1][2], t1[1]);
			}

			Fp2::sub(ret[0][2], t2[0], a[0][2]);
			Fp2::add(ret[0][2], ret[0][2], ret[0][2]);
			Fp2::add(ret[0][2], ret[0][2], t2[0]);

			Fp2::mul_by_u_plus_1(t2[1], t2[1]);
			Fp2::add(ret[1][0], t2[1], a[1][0]);
			Fp2::add(ret[1][0], ret[1][0], ret[1][0]);
			Fp2::add(ret[1][0], ret[1][0], t2[1]);

		}
	}

	BLST_HD BLST_INLINE void frobenius_map_fp6(Fp6& ret, const Fp6& a, size_t n)
	{
		const Fp2 coeffs1[] = {  /* (u + 1)^((P^n - 1) / 3) */
		  Fp2{ Fp{ 0 },
			Fp{ TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
			  TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
			  TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741) } },
		  Fp2{
			Fp{ TO_LIMB_T(0x30f1361b798a64e8), TO_LIMB_T(0xf3b8ddab7ece5a2a),
			  TO_LIMB_T(0x16a8ca3ac61577f7), TO_LIMB_T(0xc26a2ff874fd029b),
			  TO_LIMB_T(0x3636b76660701c6e), TO_LIMB_T(0x051ba4ab241b6160) } },
		  Fp2{ Fp{ 0 }, Fp{ ONE_MONT_P_LIMBS } }
		};
		const Fp coeffs2[] = {  /* (u + 1)^((2P^n - 2) / 3) */
		  {   TO_LIMB_T(0x890dc9e4867545c3), TO_LIMB_T(0x2af322533285a5d5),
			  TO_LIMB_T(0x50880866309b7e2c), TO_LIMB_T(0xa20d1b8c7e881024),
			  TO_LIMB_T(0x14e4f04fe2db9068), TO_LIMB_T(0x14e56d3f1564853a)   },
		  {   TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
			  TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
			  TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)   },
		  {   TO_LIMB_T(0x43f5fffffffcaaae), TO_LIMB_T(0x32b7fff2ed47fffd),
			  TO_LIMB_T(0x07e83a49a2e99d69), TO_LIMB_T(0xeca8f3318332bb7a),
			  TO_LIMB_T(0xef148d1ea0f4c069), TO_LIMB_T(0x040ab3263eff0206)   }
		};

		frobenius_map_fp2(ret[0], a[0], n);
		frobenius_map_fp2(ret[1], a[1], n);
		frobenius_map_fp2(ret[2], a[2], n);
		--n;    /* implied ONE_MONT_P at index 0 */
		Fp2::mul(ret[1], ret[1], coeffs1[n]);
		Fp::mul(ret[2][0], ret[2][0], coeffs2[n]);
		Fp::mul(ret[2][1], ret[2][1], coeffs2[n]);
	}

	BLST_HD FP12_CPP_INLINE void Fp12::frobenius_map(Fp12& ret, const Fp12& a, size_t n)
	{
		const Fp2 coeffs[] = {  /* (u + 1)^((P^n - 1) / 6) */
		  Fp2{
			Fp{ TO_LIMB_T(0x07089552b319d465), TO_LIMB_T(0xc6695f92b50a8313),
			  TO_LIMB_T(0x97e83cccd117228f), TO_LIMB_T(0xa35baecab2dc29ee),
			  TO_LIMB_T(0x1ce393ea5daace4d), TO_LIMB_T(0x08f2220fb0fb66eb) },
			Fp{ TO_LIMB_T(0xb2f66aad4ce5d646), TO_LIMB_T(0x5842a06bfc497cec),
			  TO_LIMB_T(0xcf4895d42599d394), TO_LIMB_T(0xc11b9cba40a8e8d0),
			  TO_LIMB_T(0x2e3813cbe5a0de89), TO_LIMB_T(0x110eefda88847faf) } },
		  Fp2{
			Fp{ TO_LIMB_T(0xecfb361b798dba3a), TO_LIMB_T(0xc100ddb891865a2c),
			  TO_LIMB_T(0x0ec08ff1232bda8e), TO_LIMB_T(0xd5c13cc6f1ca4721),
			  TO_LIMB_T(0x47222a47bf7b5c04), TO_LIMB_T(0x0110f184e51c5f59) } },
		  Fp2{
			Fp{ TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
			  TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
			  TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8) },
			Fp{ TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
			  TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
			  TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) } },
		};

		frobenius_map_fp6(ret[0], a[0], n);
		frobenius_map_fp6(ret[1], a[1], n);
		--n;    /* implied ONE_MONT_P at index 0 */
		Fp2::mul(ret[1][0], ret[1][0], coeffs[n]);
		Fp2::mul(ret[1][1], ret[1][1], coeffs[n]);
		Fp2::mul(ret[1][2], ret[1][2], coeffs[n]);
	}



	BLST_HD FP12_CPP_INLINE int Fp12::equal(const Fp12& a, const Fp12& b)
	{
		return (int)vec_is_equal(&a, &b, sizeof(Fp12));
	}

	BLST_HD FP12_CPP_INLINE int Fp12::is_one(const Fp12& a)
	{
		return
			(int)(vec_is_equal(
				&a[0][0],
				&get_BLST_Rx().p2, sizeof(a[0][0])) &
				vec_is_zero(&a[0][1], sizeof(Fp12) - sizeof(a[0][0])));
	}


	BLST_HD FP12_CPP_INLINE const Fp12& Fp12::one()
	{
		return *(Fp12*)get_BLST_Rx().p12;
	}


	BLST_HD FP12_CPP_INLINE void blst_fp12_one_copy(Fp12& ret)
	{
		char* src = (char*)get_BLST_Rx().p12;
		char* dst = (char*)&ret;
		int size = sizeof(Fp12);
		for (int i = 0; i < size; ++i)
			dst[i] = src[i];
	}

	BLST_HD FP12_CPP_INLINE void Fp12::from_bendian(Fp12& ret, const unsigned char serialized[12 * 48]) {
		// Clear the output first
		for (int i = 0; i < sizeof(Fp12); ++i)
			((unsigned char*)&ret)[i] = 0;

		// Fp12& is structured as ret[2][3][2] representing:
		// - 2 components of Fp6 (for Fp12 = Fp6[w])
		// - 3 components of Fp2 each (for Fp6 = Fp2[v])  
		// - 2 components of Fp each (for Fp2 = Fp[u])

		size_t byte_offset = 0;
		for (size_t fp6_comp = 0; fp6_comp < 3; fp6_comp++) {      // Fp6 components  
			for (size_t fp12_comp = 0; fp12_comp < 2; fp12_comp++) {       // Fp12 components
				for (size_t fp2_comp = 0; fp2_comp < 2; fp2_comp++) {  // Fp2 components
					Fp r;
					// Convert from big-endian bytes to Fp (limbs)
					limbs_from_be_bytes(r.data, &serialized[byte_offset], 48);
					Fp::to_mont(ret[fp12_comp][fp6_comp][fp2_comp], r);
					byte_offset += 48;
				}
			}
		}
	}

	BLST_HD FP12_CPP_INLINE void Fp12::to_bendian(unsigned char ret[48 * 12], const Fp12& a)
	{
		size_t i, j;
		Fp out;

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 2; j++) {
				for (int k = 0; k < 2; ++k)
				{
					Fp::from_mont(out, a[j][i][k]);
					be_bytes_from_limbs(ret, out.data, sizeof(Fp));
					ret += 48;
				}
			}
		}
	}

	BLST_HD FP12_CPP_INLINE size_t blst_fp12_sizeof(void)
	{
		return sizeof(Fp12);
	}


	BLST_HD FP12_CPP_INLINE void Fp12::exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits)
	{
		// For bits=0 or if exponent is 0, return 1
		bool_t is_zero = 1;
		for (size_t i = 0; i < (bits + 8 * sizeof(byte) - 1) / (8 * sizeof(byte)); i++) {
			if (exp[i] != 0) {
				is_zero = 0;
				break;
			}
		}

		if (bits == 0 || is_zero) {
			ret = Fp12::one();
			return;
		}

		// Make a copy of the base to prevent issues if a == ret
		Fp12 base;
		//memcpy(base, a, sizeof(Fp12));
		for (int i = 0; i < sizeof(Fp12); ++i)
			((unsigned char*)&base)[i] = ((unsigned char*)&a)[i];

		// Start with 1
		ret = Fp12::one();

		// Process bits from least to most significant
		for (size_t i = 0; i < bits; i++) {
			// If the current bit is set
			size_t limb_idx = i / (8 * sizeof(byte));
			size_t bit_idx = i % (8 * sizeof(byte));

			if ((exp[limb_idx] >> bit_idx) & 1) {
				Fp12::mul(ret, ret, base);
			}

			// Square the base for the next bit position (if not the last bit)
			if (i < bits - 1) {
				Fp12::sqr(base, base);
			}


		}
	}


	BLST_HD FP12_CPP_INLINE void Fp12::cyclotomic_exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits)
	{
		// Optimization: Dispatch to GLV method for 256-bit scalars.
		// GLV decomposition provides ~2x speedup by splitting the scalar 
		// and using the Frobenius endomorphism.
		// We restrict this to cases where bits <= 256 to ensure the scalar 
		// fits within the implementation's splitting bounds.
		if (bits > 128 && bits <= 256) {
			// GLV implementation handles variable bit lengths up to 256 internally 
			// if passed as a 256-bit container, but explicit dispatch usually 
			// implies a full-width scalar or one that fits in 256-bit words.
			limb_t exp_256[4] = { 0 };
			size_t byte_len = (bits + 7) / 8;
			memcpy(exp_256, exp, byte_len);
			//for (size_t i = 0; i < byte_len; i++) {
			//	exp_256[i / sizeof(limb_t)] |= ((limb_t)exp[i]) << ((i % sizeof(limb_t)) * 8);
			//}
			if (bits % 8 != 0) {
				// Create mask for valid bits in the top byte
				// e.g., bits % 8 == 1 -> mask 0x01
				uint8_t mask = (1 << (bits % 8)) - 1;
				// Apply mask to the specific byte in the limb array
				((uint8_t*)exp_256)[byte_len - 1] &= mask;
			}

			cyclotomic_exp_glv4(ret, a, exp_256);
		}
		else
		{
			// Fallback to basic square-and-multiply for other sizes
			cyclotomic_exp_basic(ret, a, exp, bits);
		}
	}

	BLST_HD FP12_CPP_INLINE void Fp12::cyclotomic_exp_basic(Fp12& ret, const Fp12& a, const byte* exp, size_t bits)
	{
		// 1. Find the first set bit (scanning MSB to LSB)
		int i = (int)bits - 1;
		while (i >= 0 && !((exp[i / 8] >> (i % 8)) & 1)) {
			i--;
		}

		// Handle the case where the exponent is 0
		if (i < 0) {
			ret = Fp12::one();
			return;
		}

		// 2. Initialize ret with a (the first bit)
		// This saves one mul() and prevents needing an identity '1' in registers
		ret = a;
		i--;

		// 3. Main Loop
		for (; i >= 0; i--) {
			// We use cyclotomic squaring because we are in GT
			Fp12::cyclotomic_sqr(ret, ret);

			if ((exp[i / 8] >> (i % 8)) & 1) {
				// By using 'a' directly here, the compiler can potentially 
				// load it from constant memory/L1 cache instead of 
				// keeping it in registers for the whole function.
				Fp12::mul(ret, ret, a);
			}
		}
	}


	BLST_HD FP12_CPP_INLINE void Fp12::glv_split_2(
		limb_t k1_out[4],
		limb_t k2_out[4],
		const limb_t k[4])
	{
		// Constants
		const limb_t d0 = 0x0000000100000000ULL;
		const limb_t d1 = 0xac45a4010001a402ULL;
		const limb_t mu0 = 0x63f6e522f6cfee2eULL;
		const limb_t mu1 = 0x7c6becf1e01faaddULL;

		// q_approx = (k >> 128) + ((k * mu_lo) >> 256)
		// 1. Term A = k >> 128 = { k[2], k[3], 0, 0 }
		// 2. Term B = (k * mu_lo) >> 256. 
		//    We calculate k[0..3] * mu[0..1], result 6 limbs, keep limbs [4,5].

		limb_t p[6] = { 0 };

		// k * mu0
		{
			limb_t carry = 0;
			for (int i = 0; i < 4; i++) {
				wide_limb_t w;
				cc_mul_wide_limb_limb(&w, k[i], mu0);
				cc_add_wide_wide_limb(&w, w, carry);
				p[i] = w.lo;
				carry = w.hi;
			}
			p[4] = carry;
		}

		// k * mu1, shift offset 1
		{
			limb_t carry = 0;
			for (int i = 0; i < 4; i++) {
				wide_limb_t w;
				cc_mul_wide_limb_limb(&w, k[i], mu1);
				// add to accumulator p[i+1]
				wide_limb_t sum;
				cc_add_wide_limb_limb(&sum, w.lo, p[i + 1]);
				cc_add_wide_wide_limb(&sum, sum, carry);

				p[i + 1] = sum.lo;
				// new carry is w.hi + sum.hi
				carry = w.hi + sum.hi;
			}
			p[5] = carry;
		}

		// Sum q = Term A + Term B
		// q[0] = k[2] + p[4], q[1] = k[3] + p[5]
		limb_t q[4] = { 0 };

		{
			wide_limb_t w;
			cc_add_wide_limb_limb(&w, k[2], p[4]);
			q[0] = w.lo;

			wide_limb_t w2;
			cc_add_wide_limb_limb(&w2, k[3], p[5]);
			cc_add_wide_wide_limb(&w2, w2, w.hi); // carry
			q[1] = w2.lo;
			q[2] = w2.hi; // q might overflow 128 bits slightly
		}

		// Calculate Remainder r = k - q*d
		// We only need low 256 bits of q*d.
		limb_t prod[4] = { 0 };

		auto mul_add_d = [&](limb_t q_limb, int offset) {
			if (q_limb == 0) return;
			wide_limb_t w;
			limb_t carry = 0;

			// d0
			cc_mul_wide_limb_limb(&w, d0, q_limb);
			cc_add_wide_wide_limb(&w, w, prod[offset]); // add to existing
			prod[offset] = w.lo;
			carry = w.hi;

			// d1
			if (offset + 1 < 4) {
				wide_limb_t w2;
				cc_mul_wide_limb_limb(&w2, d1, q_limb);
				cc_add_wide_wide_limb(&w2, w2, prod[offset + 1]);
				cc_add_wide_wide_limb(&w2, w2, carry);
				prod[offset + 1] = w2.lo;
				carry = w2.hi;
			}

			// propagate
			for (int i = offset + 2; i < 4 && carry; i++) {
				wide_limb_t s;
				cc_add_wide_limb_limb(&s, prod[i], carry);
				prod[i] = s.lo;
				carry = s.hi;
			}
			};

		mul_add_d(q[0], 0);
		mul_add_d(q[1], 1);
		mul_add_d(q[2], 2);

		// r = k - prod
		limb_t r[4];
		limb_t borrow = 0;
		for (int i = 0; i < 4; i++) {
			limb_t diff = k[i] - prod[i] - borrow;
			borrow = (k[i] < prod[i]) | ((k[i] == prod[i]) & borrow);
			r[i] = diff;
		}

		// Correction loop
		// Because q is floor estimate (or slightly less), q*d <= k is expected.
		// So r should be positive (borrow == 0).
		// If borrow == 1, then prod > k, meaning q was too BIG.
		// This happens if our estimate was slightly over?
		// Actually q_approx <= k/d generally.
		// But let's handle both cases.

		limb_t d_full[4] = { d0, d1, 0, 0 };

		auto ge = [&](const limb_t A[4], const limb_t B[4]) {
			for (int i = 3; i >= 0; i--) if (A[i] != B[i]) return A[i] > B[i];
			return true;
			};

		if (borrow) {
			// Wraparound means negative result (k < prod). q too big.
			// r currently holds (k - prod) mod 2^256, which is large.
			// We need to add d back to r and decrement q until r is small positive.
			// Effectively: while (r is "negative") { r += d; q--; }
			while (borrow) {
				// r += d
				limb_t carry = 0;
				for (int i = 0; i < 4; i++) {
					wide_limb_t w;
					cc_add_wide_limb_limb(&w, r[i], d_full[i]);
					cc_add_wide_wide_limb(&w, w, carry);
					r[i] = w.lo;
					carry = w.hi;
				}
				// If carry out is 1, we wrapped back to positive range
				if (carry) borrow = 0;

				// q--
				limb_t sub_b = 1;
				for (int i = 0; i < 4; i++) {
					limb_t res = q[i] - sub_b;
					sub_b = (q[i] < sub_b);
					q[i] = res;
				}
			}
		}
		else {
			// r is positive. Check if r >= d (q too small).
			while (ge(r, d_full)) {
				// r -= d
				limb_t b_inner = 0;
				for (int i = 0; i < 4; i++) {
					limb_t diff = r[i] - d_full[i] - b_inner;
					b_inner = (r[i] < d_full[i]) | ((r[i] == d_full[i]) & b_inner);
					r[i] = diff;
				}

				// q++
				limb_t c_inner = 1;
				for (int i = 0; i < 4; i++) {
					wide_limb_t w;
					cc_add_wide_limb_limb(&w, q[i], c_inner);
					q[i] = w.lo;
					c_inner = w.hi;
				}
			}
		}

		// Copy out
		for (int i = 0; i < 4; i++) {
			k1_out[i] = r[i];
			k2_out[i] = q[i];
		}
	}


	// Overload for pre-split scalars
	BLST_HD FP12_CPP_INLINE void Fp12::cyclotomic_exp_glv(Fp12& out, const Fp12& base, const limb_t* k1, const limb_t* k2)
	{
		// Precompute window table
		// T[0] = 1 (implied)
		// T[1] = g^1 * phi(g)^0 = g
		// T[2] = g^0 * phi(g)^1 = phi^2(g)
		// T[3] = g^1 * phi(g)^1 = g * phi^2(g)

		Fp12 g_phi2;
		Fp12::frobenius_map(g_phi2, base, 2);

		// We need T[3] = base * g_phi2
		Fp12 T3;
		Fp12::mul(T3, base, g_phi2);

		// Initialize accumulator with identity? 
		// Or better: find MSB and start there to save squarings.
		// max bits of k1, k2 is approx 128.

		vec_copy(&out, &Fp12::one(), sizeof(Fp12));

		// Scan bits from 128 down to 0
		// We can just run fixed 129 iterations for constant-time-ishness
		// (though split logic had data dependent loops, so full CT requires more work)

		for (int i = 128; i >= 0; --i) {
			Fp12::cyclotomic_sqr(out, out);

			unsigned b1 = (unsigned)((k1[i / 64] >> (i % 64)) & 1);
			unsigned b2 = (unsigned)((k2[i / 64] >> (i % 64)) & 1);
			unsigned idx = b1 | (b2 << 1); // 0..3

			if (idx == 1) Fp12::mul(out, out, base);
			else if (idx == 2) Fp12::mul(out, out, g_phi2);
			else if (idx == 3) Fp12::mul(out, out, T3);
		}
	}

	// Original GLV exponentiation calling the splitter
	BLST_HD FP12_CPP_INLINE void Fp12::cyclotomic_exp_glv(
		Fp12& out, const Fp12& base, const byte* scalar_bytes)
	{
		// Convert input bytes to limbs
		limb_t k[4];
		// Assume LE bytes
		// Simple load if aligned, else cautious load
		const u64* k64 = (const u64*)scalar_bytes; // Warning: alignment
		// Safe copy:
		for (int i = 0; i < 4; i++) {
			// primitive manual load for portability/safety
			u64 w = 0;
			for (int b = 0; b < 8; b++) w |= ((u64)scalar_bytes[i * 8 + b]) << (b * 8);
			k[i] = (limb_t)w;
		}

		limb_t k1[4], k2[4];
		glv_split_2(k1, k2, k);

		// Delegate to core logic
		cyclotomic_exp_glv(out, base, k1, k2);
	}

	// Helper for safe bit access
	BLST_HD inline unsigned get_scalar_bit(const limb_t* k, int i) {
		if (i < 0) return 0;
		// Assumes 64-bit limbs.
		return (unsigned)((k[i / 64] >> (i % 64)) & 1);
	}

	BLST_HD GLV_JSF_Exp Fp12::glv_recode_jsf(const limb_t* k1, const limb_t* k2)
	{
		GLV_JSF_Exp result = { 0 }; // Zero-initialize for safe bitwise OR in set()
		unsigned d1 = 0, d2 = 0; // Carries for k1, k2
		unsigned u, v;
		int i = 0;
		result.sched_len = 0;

		// Map scalars to signed digits {-1, 0, 1}
		// Loop limit 132 covers 129-bit input + worst-case carry propagation
		while (i < 133) {
			unsigned k1_bit = get_scalar_bit(k1, i);
			unsigned k2_bit = get_scalar_bit(k2, i);

			// Lookahead bits for JSF decision
			unsigned k1_next = get_scalar_bit(k1, i + 1);
			unsigned k2_next = get_scalar_bit(k2, i + 1);

			u = k1_bit + d1;
			v = k2_bit + d2;

			int j1, j2;

			// JSF Logic for row 1
			if ((u & 1) == 0) {
				j1 = 0;
			}
			else {
				// Cohen's JSF Logic:
				// If u == 3 (mod 4) OR (u == 1 (mod 4) AND v == 2 (mod 4)), choose -1
				// Note: v % 4 checks lookahead of v.
				unsigned u_mod_4 = (u + (k1_next << 1)) & 3;
				unsigned v_mod_4 = (v + (k2_next << 1)) & 3;

				if (u_mod_4 == 3 || (u_mod_4 == 1 && v_mod_4 == 2)) j1 = -1;
				else j1 = 1;
			}

			// JSF Logic for row 2
			if ((v & 1) == 0) {
				j2 = 0;
			}
			else {
				unsigned u_mod_4 = (u + (k1_next << 1)) & 3;
				unsigned v_mod_4 = (v + (k2_next << 1)) & 3;

				if (v_mod_4 == 3 || (v_mod_4 == 1 && u_mod_4 == 2)) j2 = -1;
				else j2 = 1;
			}

			// Bounds check for safety before writing
			if (result.sched_len < 132) {
				result.set(result.sched_len, j1, j2);
				result.sched_len++;
			}

			// Calculate next carries: d = (u - j) / 2
			d1 = (u - j1) >> 1;
			d2 = (v - j2) >> 1;

			i++;

			// Quit early if we are past bit 130 and have no carries left
			if (i >= 130 && d1 == 0 && d2 == 0) break;
		}

		// Trim leading zeros to find real MSB
		int len = result.sched_len;
		while (len > 0 && result.get(len - 1).d1 == 0 && result.get(len - 1).d2 == 0) {
			len--;
		}
		result.sched_len = (uint8_t)len;

		return result;
	}

	//// Full GLV with JSF recoding
	//// Phase 2: Lean Execution Loop (Split API)
	//BLST_HD void Fp12::cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, int8_t jsf[132][2], int len)
	//{
	//	if (len == 0) {
	//		vec_copy(&out, &Fp12::one(), sizeof(Fp12));
	//		return;
	//	}

	//	// Precompute lookup table
	//	// 0:unused, 1:P1, 2:P2, 3:P1*P2, 4:P1*conj(P2)
	//	Fp12 P2, P1_P2, P1_nP2;
	//	Fp12::frobenius_map(P2, base, 2);
	//	Fp12::mul(P1_P2, base, P2);
	//	{
	//		Fp12 P2_conj = P2;
	//		Fp12::conjugate(P2_conj);
	//		Fp12::mul(P1_nP2, base, P2_conj);
	//	}

	//	const Fp12* LUT[5] = { nullptr, &base, &P2, &P1_P2, &P1_nP2 };

	//	// Initialize with MSB
	//	len--; // Point to last valid index
	//	{
	//		int u = jsf[len][0];
	//		int v = jsf[len][1];
	//		int idx = 0;
	//		bool invert = false;

	//		if (u == 0) {
	//			idx = 2; // (0, +/-1) -> P2
	//			if (v < 0) invert = true;
	//		}
	//		else if (v == 0) {
	//			idx = 1; // (+/-1, 0) -> P1
	//			if (u < 0) invert = true;
	//		}
	//		else if (u == v) {
	//			idx = 3; // (1, 1) or (-1, -1) -> P1*P2
	//			if (u < 0) invert = true;
	//		}
	//		else {
	//			idx = 4; // (1, -1) or (-1, 1) -> P1*conj(P2)
	//			if (u < 0) invert = true;
	//		}

	//		vec_copy(&out, LUT[idx], sizeof(Fp12));
	//		if (invert) Fp12::conjugate(out);
	//	}

	//	// Main Loop
	//	for (len--; len >= 0; --len) {
	//		Fp12::cyclotomic_sqr(out, out);

	//		int u = jsf[len][0];
	//		int v = jsf[len][1];

	//		if (u == 0 && v == 0) continue;

	//		int idx = 0;
	//		bool invert = false;

	//		if (u == 0) {
	//			idx = 2;
	//			if (v < 0) invert = true;
	//		}
	//		else if (v == 0) {
	//			idx = 1;
	//			if (u < 0) invert = true;
	//		}
	//		else if (u == v) {
	//			idx = 3;
	//			if (u < 0) invert = true;
	//		}
	//		else {
	//			idx = 4;
	//			if (u < 0) invert = true;
	//		}

	//		if (invert) {
	//			Fp12 tmp = *LUT[idx];
	//			Fp12::conjugate(tmp);
	//			Fp12::mul(out, out, tmp);
	//		}
	//		else {
	//			Fp12::mul(out, out, *LUT[idx]);
	//		}
	//	}
	//}

	// Optimized for GPU/Stack-heavy environments: 
	// Recodes on the fly to avoid memory bandwidth of passing JSF array.
	BLST_HD void Fp12::cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const limb_t* k1, const limb_t* k2)
	{
		cyclotomic_exp_glv_jsf(out, base, glv_recode_jsf(k1, k2));
	}

	// Struct for pair return to avoid std::pair dependency on device
	struct OpSelect { const Fp12* ptr; bool inv; };
	
	BLST_HD auto get_op(int u, int v, const Fp12& base, const Fp12& P2, const Fp12& P1_P2, const Fp12& P1_nP2) {
		const Fp12* op = nullptr;
		bool inv = false;
		if (u == 0) { op = &P2; if (v < 0) inv = true; }
		else if (v == 0) { op = &base; if (u < 0) inv = true; }
		else if (u == v) { op = &P1_P2; if (u < 0) inv = true; }
		else { op = &P1_nP2; if (u < 0) inv = true; }
		return OpSelect{ op, inv };
	};

	BLST_HD void Fp12::cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const GLV_JSF_Exp& jsf)
	{
		int len = jsf.sched_len;
		if (len == 0) {
			vec_copy(&out, &Fp12::one(), sizeof(Fp12));
			return;
		}

		// 2. Precompute Lookup Points (Stack)
		// We compute mixed points to reduce the loop to 1 add per step max.
		Fp12 P2, P1_P2, P1_nP2;

		// P2 = phi^2(P1)
		Fp12::frobenius_map(P2, base, 2);

		// P1_P2 = P1 * P2
		Fp12::mul(P1_P2, base, P2);

		// P1_nP2 = P1 * conj(P2)
		Fp12::mul(P1_nP2, base, P2, true);

		// 3. Init with MSB
		len--;

		// Initialization
		{
			auto pair = jsf.get(len);
			OpSelect sel = get_op(pair.d1, pair.d2, base, P2, P1_P2, P1_nP2);

			vec_copy(&out, sel.ptr, sizeof(Fp12));
			if (sel.inv) Fp12::conjugate(out);
		}

		// 4. Main Loop
		for (len--; len >= 0; --len) {
			Fp12::cyclotomic_sqr(out, out);

			auto pair = jsf.get(len);
			int u = pair.d1;
			int v = pair.d2;

			if (u == 0 && v == 0)
				continue;
			
			// Minimized divergence: 
			// No stack allocations, no "conjugate" function call.
			// Just a boolean flag passed to mul.
			OpSelect sel = get_op(u, v, base, P2, P1_P2, P1_nP2);
			Fp12::mul(out, out, *sel.ptr, sel.inv);
		}
	}

	BLST_HD FP12_CPP_INLINE bool_t Fp12::is_cyclotomic(const Fp12& f)
	{
		Fp12 a, b;

		Fp12::frobenius_map(a, f, 2);
		Fp12::frobenius_map(b, a, 2);
		Fp12::mul(b, b, f);

		return vec_is_equal(&a, &b, sizeof(a));
	}


	BLST_HD FP12_CPP_INLINE void Fp12::mul_n_sqr(Fp12& ret, const Fp12& a, size_t n)
	{
		Fp12::mul(ret, ret, a);
		while (n--)
			Fp12::cyclotomic_sqr(ret, ret);
	}

	BLST_HD FP12_CPP_INLINE void Fp12::raise_to_z_div_by_2(Fp12& ret, const Fp12& a)
	{
		Fp12::cyclotomic_sqr(ret, a);                /* 0x2                  */
		Fp12::mul_n_sqr(ret, a, 2);                       /* ..0xc                */
		Fp12::mul_n_sqr(ret, a, 3);                       /* ..0x68               */
		Fp12::mul_n_sqr(ret, a, 9);                       /* ..0xd200             */
		Fp12::mul_n_sqr(ret, a, 32);                      /* ..0xd20100000000     */
		Fp12::mul_n_sqr(ret, a, 16 - 1);                    /* ..0x6900800000008000 */
		Fp12::conjugate(ret);                /* account for z being negative */
	}

	BLST_HD FP12_CPP_INLINE void Fp12::raise_to_z(Fp12& a, const Fp12& b)
	{
		raise_to_z_div_by_2(a, b);
		Fp12::cyclotomic_sqr(a, a);
	}

	BLST_HD FP12_CPP_INLINE int Fp12::in_group(const Fp12& f)
	{
		Fp12 a, b;

		if (vec_is_zero(&f, sizeof(Fp12)) || !is_cyclotomic(f))
			return 0;

		Fp12::frobenius_map(a, f, 1);
		raise_to_z(b, f);

		return (int)vec_is_equal(&a, &b, sizeof(a));
	}



	/*
	 * Adaptation from <zkcrypto>/pairing/src/bls12_381/mod.rs
	 */
	BLST_HD FP12_CPP_INLINE void Fp12::final_exp(Fp12& ret, const Fp12& f)
	{
		Fp12 y0, y1, y2, y3;

		vec_copy(&y1, &f, sizeof(y1));
		Fp12::conjugate(y1);
		Fp12::inverse(y2, f);
		Fp12::mul(ret, y1, y2);
		Fp12::frobenius_map(y2, ret, 2);
		Fp12::mul(ret, ret, y2);

		Fp12::cyclotomic_sqr(y0, ret);
		Fp12::raise_to_z(y1, y0);
		Fp12::raise_to_z_div_by_2(y2, y1);
		vec_copy(&y3, &ret, sizeof(y3));
		Fp12::conjugate(y3);
		Fp12::mul(y1, y1, y3);
		Fp12::conjugate(y1);
		Fp12::mul(y1, y1, y2);
		Fp12::raise_to_z(y2, y1);
		Fp12::raise_to_z(y3, y2);
		Fp12::conjugate(y1);
		Fp12::mul(y3, y3, y1);
		Fp12::conjugate(y1);
		Fp12::frobenius_map(y1, y1, 3);
		Fp12::frobenius_map(y2, y2, 2);
		Fp12::mul(y1, y1, y2);
		Fp12::raise_to_z(y2, y3);
		Fp12::mul(y2, y2, y0);
		Fp12::mul(y2, y2, ret);
		Fp12::mul(y1, y1, y2);
		Fp12::frobenius_map(y2, y3, 1);
		Fp12::mul(ret, y1, y2);
	}


	// 128-bit / 64-bit division: (hi, lo) / div
	// Returns quotient, puts remainder in *rem
	BLST_HD FP12_CPP_INLINE static uint64_t udiv128_64(uint64_t hi, uint64_t lo, uint64_t div, uint64_t* rem) {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__CUDA_ARCH__)
		uint64_t r;
		uint64_t q = _udiv128(hi, lo, div, &r);
		*rem = r;
		return q;
#elif defined(__x86_64__) || defined(__aarch64__)
		unsigned __int128 p = ((unsigned __int128)hi << 64) | lo;
		*rem = (uint64_t)(p % div);
		return (uint64_t)(p / div);
#else
		// Fallback for purely portable arithmetic (bit-by-bit or half-word)
		// Simplified bitwise restoring division for test purposes if intrinsics missing
		// Slow but functional.
		if (hi >= div) { /* overflow case, strict 128/64 implies q fits in 64, so hi < div */ }
		uint64_t r = hi;
		uint64_t q = 0;
		for (int i = 63; i >= 0; i--) {
			r = (r << 1) | ((lo >> i) & 1);
			if (r >= div) {
				r -= div;
				q |= (1ULL << i);
			}
		}
		*rem = r;
		return q;
#endif
	}

	// Helper: 256-bit by 64-bit division/remainder.
	// Computes q = n / d, r = n % d.
	BLST_HD FP12_CPP_INLINE static void div_rem_u256_u64(limb_t q[4], uint64_t& r, const limb_t n[4], uint64_t d)
	{
		uint64_t rem = 0;
		for (int i = 3; i >= 0; i--) {
			q[i] = udiv128_64(rem, n[i], d, &rem);
		}
		r = rem;
	}
	// BLS12-381 Parameter z (absolute value)
	static const uint64_t BLS_Z = 0xd201000000010000ULL;
	// Split 256-bit scalar k into 4 64-bit signed scalars k0..k3
	// such that k = k0 + k1*x + k2*x^2 + k3*x^3 mod r
	// where x is the BLS parameter (negative).
	BLST_HD FP12_CPP_INLINE
		void Fp12::gt_split_scalar_4d(int64_t k_out[4], const limb_t k_in[4])
	{
		limb_t temp_k[4];
		vec_copy(temp_k, k_in, sizeof(limb_t) * 4);

		// Base-z expansion
		for (int i = 0; i < 4; i++) {
			uint64_t r = 0;
			limb_t q[4];

			// If temp_k is 0, we are done
			if ((temp_k[0] | temp_k[1] | temp_k[2] | temp_k[3]) == 0) {
				k_out[i] = 0;
				continue;
			}

			div_rem_u256_u64(q, r, temp_k, BLS_Z);

			// Check balance
			bool carry = false;
			if (r > (BLS_Z / 2)) {
				carry = true;
				k_out[i] = (int64_t)r - (int64_t)BLS_Z;
			}
			else {
				k_out[i] = (int64_t)r;
			}

			// Update temp_k with q
			// If carry, we need to add 1 to q. 
			// Check for overflow in q (very rare given q is roughly 256-64 bits each step)
			if (carry) {
				int j = 0;
				while (j < 4) {
					q[j]++;
					if (q[j] != 0) break;
					j++;
				}
			}
			vec_copy(temp_k, q, sizeof(limb_t) * 4);
		}

		// Handle carry-out (k4) which corresponds to z^4 = u^4 = (r + u^2 - 1)
		// Since we are mod r, this is u^2 - 1.
		// Effect on k: add k4*u^2 - k4*1
		// Base-z decomposition: k = k0 + k1(-z) + k2(-z)^2 + k3(-z)^3 + k4(-z)^4
		// k = k0 - k1*z + k2*z^2 - k3*z^3 + k4*z^4
		// Coeff of z^2 (u^2) is increased by k4.
		// Coeff of z^0 (1) is decreased by k4.
		// temp_k is the remaining quotient, used as k4.
		if (temp_k[0] != 0) {
			k_out[2] += (int64_t)temp_k[0];
			k_out[0] -= (int64_t)temp_k[0];
		}

		// Adjust for x = -z
		k_out[1] = -k_out[1];
		k_out[3] = -k_out[3];
	}

	BLST_HD FP12_CPP_INLINE
	void Fp12::cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const limb_t* scalar)
	{


#ifndef NDEBUG
		// Debug assertion: Scalar must be strictly less than group order r.
		// GLV decomposition optimized for curve operations assumes valid Fr scalars.
		{
			const limb_t* r = get_BLST_r().data;
			int cmp = 0;
			for (int i = 3; i >= 0; i--) {
				if (scalar[i] > r[i]) { cmp = 1; break; }
				if (scalar[i] < r[i]) { cmp = -1; break; }
			}
			// If equal (0) or greater (1), assert fails
			BLST_ASSERT(cmp < 0 && "Scalar must be reduced ( < r ) for GLV optimization");
		}
#endif

		// 1. Decompose Scalar
		int64_t k[4];
		gt_split_scalar_4d(k, scalar);

		// 2. Prepare Basis Points {P0, P1, P2, P3}
		uint64_t abs_k[4];
		Fp12 T[16];

		T[1] = base; // P0
		Fp12::frobenius_map(T[2], base, 1); // P1
		Fp12::frobenius_map(T[4], base, 2); // P2
		Fp12::frobenius_map(T[8], base, 3); // P3

		// Handle signs
		for (int i = 0; i < 4; ++i) {
			if (k[i] < 0) {
				abs_k[i] = (uint64_t)(-k[i]);
				Fp12::conjugate(T[1 << i]);
			}
			else {
				abs_k[i] = (uint64_t)k[i];
			}
		}

		// 3. Precompute Lookup Table (16 points)
		// 0000 -> 1
		vec_copy(&T[0], &Fp12::one(), sizeof(Fp12));


		// Fill rest
		Fp12::mul(T[3], T[1], T[2]);  // 1|2
		Fp12::mul(T[5], T[1], T[4]);  // 1|4
		Fp12::mul(T[6], T[2], T[4]);  // 2|4
		Fp12::mul(T[7], T[3], T[4]);  // 1|2|4

		Fp12::mul(T[9], T[1], T[8]);  // 1|8
		Fp12::mul(T[10], T[2], T[8]); // 2|8
		Fp12::mul(T[11], T[3], T[8]); // 1|2|8
		Fp12::mul(T[12], T[4], T[8]); // 4|8
		Fp12::mul(T[13], T[5], T[8]); // 1|4|8
		Fp12::mul(T[14], T[6], T[8]); // 2|4|8
		Fp12::mul(T[15], T[7], T[8]); // 1|2|4|8

		// 4. Main 64-bit Loop
		// Note limits: max k_i is roughly z/2 + carry ~ 64 bits.

		int msb = 63;
		// Find first non-zero logic could be added here, 
		// but 64 constant iterations is fine.

		// Initialize accumulator ??
		// Standard MSM: Start with identity, or find first set bit.
		// Let's find real MSB to avoid initial identity multiplies.
		uint64_t combined = abs_k[0] | abs_k[1] | abs_k[2] | abs_k[3];
		while (msb >= 0 && ((combined >> msb) & 1) == 0) msb--;

		if (msb < 0) {
			vec_copy(&out, &Fp12::one(), sizeof(Fp12));
			return;
		}

		// Init with first window
		unsigned idx = 0;
		for (int i = 0; i < 4; ++i) if ((abs_k[i] >> msb) & 1) idx |= (1 << i);
		out = T[idx];

		for (int j = msb - 1; j >= 0; --j) {
			Fp12::cyclotomic_sqr(out, out);

			idx = 0;
			for (int i = 0; i < 4; ++i) if ((abs_k[i] >> j) & 1) idx |= (1 << i);

			if (idx != 0) Fp12::mul(out, out, T[idx]);
		}
	}

	BLST_HD FP12_CPP_INLINE
	GLV4_Exp Fp12::glv4_compile(const limb_t* k_in) {
		GLV4_Exp props;
		props.sign_mask = 0;
		props.sched_len = 0;

		// ... (splitting & sign logic same as before) ...
		int64_t k[4];
		gt_split_scalar_4d(k, k_in);

		uint64_t abs_k[4];
		for (int i = 0; i < 4; i++) {
			if (k[i] < 0) {
				props.sign_mask |= (1 << i);
				abs_k[i] = (uint64_t)(-k[i]);
			}
			else {
				abs_k[i] = (uint64_t)k[i];
			}
		}

		// Helper to push
		auto push = [&](uint8_t val) {
			if (props.sched_len < 72) {
				props.schedule[props.sched_len++] = val;
			}
			};

		uint64_t combined = abs_k[0] | abs_k[1] | abs_k[2] | abs_k[3];
		int msb = 63;
		while (msb >= 0 && ((combined >> msb) & 1) == 0) msb--;

		if (msb < 0) return props;

		unsigned idx = 0;
		if ((abs_k[0] >> msb) & 1) idx |= 1;
		if ((abs_k[1] >> msb) & 1) idx |= 2;
		if ((abs_k[2] >> msb) & 1) idx |= 4;
		if ((abs_k[3] >> msb) & 1) idx |= 8;

		push((uint8_t)idx);

		int sqr_run = 1;

		for (int bit = msb - 1; bit >= 0; --bit) {
			idx = 0;
			if ((abs_k[0] >> bit) & 1) idx |= 1;
			if ((abs_k[1] >> bit) & 1) idx |= 2;
			if ((abs_k[2] >> bit) & 1) idx |= 4;
			if ((abs_k[3] >> bit) & 1) idx |= 8;

			if (idx == 0) {
				sqr_run++;
			}
			else {
				while (sqr_run > 15) {
					push((15 << 4) | 0);
					sqr_run -= 15;
				}
				push((uint8_t)((sqr_run << 4) | idx));
				sqr_run = 1;
			}
		}

		sqr_run--;
		while (sqr_run > 0) {
			int chunk = (sqr_run > 15) ? 15 : sqr_run;
			push((chunk << 4) | 0);
			sqr_run -= chunk;
		}

		return props;
	}

	BLST_HD FP12_CPP_INLINE
	void Fp12::cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const GLV4_Exp& props)
	{
		// ... (Setup: Signs, Basis, Table same as before) ...
		Fp12 basis[4];
		basis[0] = base;
		Fp12 phi;
		Fp12::frobenius_map(phi, base, 1);
		basis[1] = phi;
		Fp12::frobenius_map(phi, phi, 1);
		basis[2] = phi;
		Fp12::frobenius_map(phi, phi, 1);
		basis[3] = phi;

		if (props.sign_mask & 1) Fp12::conjugate(basis[0]);
		if (props.sign_mask & 2) Fp12::conjugate(basis[1]);
		if (props.sign_mask & 4) Fp12::conjugate(basis[2]);
		if (props.sign_mask & 8) Fp12::conjugate(basis[3]);

		Fp12 T[16];
		T[0] = basis[0]; // dummy/unused
		T[1] = basis[0];
		T[2] = basis[1];
		T[4] = basis[2];
		T[8] = basis[3];

		Fp12::mul(T[3], T[1], T[2]);
		Fp12::mul(T[5], T[1], T[4]);
		Fp12::mul(T[6], T[2], T[4]);
		Fp12::mul(T[7], T[3], T[4]);

		Fp12::mul(T[9], T[1], T[8]);
		Fp12::mul(T[10], T[2], T[8]);
		Fp12::mul(T[11], T[3], T[8]);
		Fp12::mul(T[12], T[4], T[8]);
		Fp12::mul(T[13], T[5], T[8]);
		Fp12::mul(T[14], T[6], T[8]);
		Fp12::mul(T[15], T[7], T[8]);

		// Execution from fixed array
		if (props.sched_len == 0) {
			vec_copy(&out, &Fp12::one(), sizeof(Fp12));
			return;
		}

		uint8_t op = props.schedule[0];
		int idx = op & 0xF;
		out = T[idx];

		for (int i = 1; i < props.sched_len; ++i) {
			op = props.schedule[i];
			int sqr_cnt = op >> 4;
			int m_idx = op & 0xF;

			for (int s = 0; s < sqr_cnt; ++s) Fp12::cyclotomic_sqr(out, out);

			if (m_idx != 0) {
				Fp12::mul(out, out, T[m_idx]);
			}
		}
	}

}

#endif // !BLST_FP12_TOWER_H


