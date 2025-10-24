/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/curve/point.h"

namespace blst
{

 /*
  * This implementation uses explicit addition formula:
  *
  * λ = (Y₂-Y₁)/(X₂-X₁)
  * X₃ = λ²-(X₁+X₂)
  * Y₃ = λ⋅(X₁-X₃)-Y₁
  *
  * But since we don't know if we'll have to add point to itself, we need
  * to eventually resort to corresponding doubling formula:
  *
  * λ = 3X₁²/2Y₁
  * X₃ = λ²-2X₁
  * Y₃ = λ⋅(X₁-X₃)-Y₁
  *
  * The formulae use prohibitively expensive inversion, but whenever we
  * have a lot of affine points to accumulate, we can amortize the cost
  * by applying Montgomery's batch inversion approach. As a result,
  * asymptotic[!] per-point cost for addition is as small as 5M+1S. For
  * comparison, ptype##_dadd takes 8M+5S. In practice, all things
  * considered, the improvement coefficient varies from 60% to 85%
  * depending on platform and curve.
  *
  * THIS IMPLEMENTATION IS *NOT* CONSTANT-TIME. [But if there is an
  * application that requires constant time-ness, speak up!]
  */

  /*
   * Calculate λ's numerator and denominator.
   *
   * input:	A	x1	y1	-
   *		B	x2	y2	-
   * output:
   * if A!=B:	A	x1	y1	(x2-x1)*mul_acc
   *		B	x2+x1	y2-y1	(x2-x1)
   *
   * if A==B:	A	x	y	2y*mul_acc
   *		B	2x	3*x^2	2y
   *
   * if A==-B:	A	0	0	1*mul_acc
   *		B	0	3*x^2	0
   */
//#define HEAD(ptype, bits, field, field_ptr, one) \
//void ptype##_head(ptype AB[2], const field_ptr* mul_acc) \
//{ \
//    ptype *A = AB, *B = AB+1; \
//    limb inf = vec_is_zero(A, sizeof(ptype##_affine)) | \
//                 vec_is_zero(B, sizeof(ptype##_affine));  \
//    const vec##bits zero = { 0 }; \
//\
//    sub_##field(B->Z, B->X, A->X);		/* X2-X1  */ \
//    add_##field(B->X, B->X, A->X);		/* X2+X1  */ \
//    add_##field(A->Z, B->Y, A->Y);		/* Y2+Y1  */ \
//    sub_##field(B->Y, B->Y, A->Y);		/* Y2-Y1  */ \
//    if (vec_is_zero(B->Z, sizeof(B->Z))) {	/* X2==X1 */ \
//        inf = vec_is_zero(A->Z, sizeof(A->Z));	\
//        vec_select(B->X, A->Z, B->X, sizeof(B->X), inf); \
//        sqr_##field(B->Y, A->X);		\
//        mul_by_3_##field(B->Y, B->Y);		/* 3*X1^2 */ \
//        vec_copy(B->Z, A->Z, sizeof(B->Z));	/* 2*Y1   */ \
//    }						/* B->Y is numenator    */ \
//						/* B->Z is denominator  */ \
//    vec_select(A->X, B->X, A->X, sizeof(A->X), inf); \
//    vec_select(A->Y, A->Z, A->Y, sizeof(A->Y), inf); \
//    vec_select(A->Z, one,  B->Z, sizeof(A->Z), inf); \
//    vec_select(B->Z, zero, B->Z, sizeof(B->Z), inf); \
//    if (mul_acc != NULL) \
//        mul_##field(A->Z, A->Z, mul_acc);	/* chain multiplication */\
//}
//
//   /*
//	* Calculate λ and resulting coordinates.
//	*
//	* input:	A		x1			y1		-
//	*		B		x2+x1			nominator	-
//	* 		lambda		1/denominator
//	* output:	D		x3=(nom/den)^2-(x2+x1)	y3=(nom/den)(x1-x3)-y1
//	*/
//#define TAIL(ptype, bits, field, one) \
//void ptype##_tail(ptype *D, ptype AB[2], vec##bits lambda) \
//{ \
//    ptype *A = AB, *B = AB+1; \
//    vec##bits llambda; \
//    limb inf = vec_is_zero(B->Z, sizeof(B->Z)); \
//\
//    mul_##field(lambda, lambda, B->Y);		/* λ = (Y2-Y1)/(X2-X1)  */ \
//						/* alt. 3*X1^2/2*Y1     */ \
//    sqr_##field(llambda, lambda); \
//    sub_##field(D->X, llambda, B->X);		/* X3 = λ^2-X1-X2       */ \
//\
//    sub_##field(D->Y, A->X, D->X);   \
//    mul_##field(D->Y, D->Y, lambda); \
//    sub_##field(D->Y, D->Y, A->Y);		/* Y3 = λ*(X1-X3)-Y1    */ \
//\
//    vec_select(D->X, A->X, D->X, 2*sizeof(D->X), inf); \
//    vec_select(B->Z, one, B->Z, sizeof(B->Z), inf); \
//}
//
//	/*
//	 * |points[]| is volatile buffer with |X|s and |Y|s initially holding
//	 * input affine coordinates, and with |Z|s being used as additional
//	 * temporary storage [unrelated to Jacobian coordinates]. |sum| is
//	 * in-/output, initialize to infinity accordingly.
//	 */
//#define ADDITION_BTREE(prefix, ptype, bits, field, field_ptr, one) \
//HEAD(ptype, bits, field, field_ptr, one) \
//TAIL(ptype, bits, field, one) \
//void ptype##s_accumulate(ptype *sum, ptype points[], size_t n) \
//{ \
//    ptype *dst; \
//    void *mul_acc; \
//    size_t i; \
//\
//    while (n >= 16) { \
//        if (n & 1) \
//            ptype##_dadd(sum, sum, (const ptype##_affine *)points++); \
//        n /= 2; \
//        for (mul_acc = NULL, i = n; i--; mul_acc = points->Z, points += 2) \
//            ptype##_head(points, (const field_ptr*)mul_acc); \
//\
//        reciprocal_##field(points[-2].Z, points[-2].Z); /* 1/∏ Zi */ \
//\
//        for (dst = points, i = n; --i;) { \
//            dst--; points -= 2; \
//            mul_##field(points[-2].Z, points[0].Z, points[-2].Z); \
//            ptype##_tail(dst, points, points[-2].Z); \
//            mul_##field(points[-2].Z, points[0].Z, points[1].Z); \
//        } \
//        dst--; points -= 2; \
//        ptype##_tail(dst, points, points[0].Z); \
//        points = dst; \
//    } \
//    while (n--) \
//        ptype##_dadd(sum, sum, (const ptype##_affine *)points++); \
//} \
//\
//void prefix##s_add(ptype *sum, const ptype##_affine *const points[], \
//                               size_t npoints) \
//{ \
//    const size_t stride = SCRATCH_LIMIT / sizeof(ptype); \
//    ptype *scratch = (ptype*)alloca((npoints > stride ? stride : npoints) * \
//                            sizeof(ptype)); \
//    const ptype##_affine *point = NULL; \
//\
//    vec_zero(sum, sizeof(*sum)); \
//    while (npoints) { \
//        size_t i, j = npoints > stride ? stride : npoints; \
//        for (i=0; i<j; i++) { \
//            point = *points ? *points++ : point+1; \
//            vec_copy(&scratch[i], point, sizeof(*point)); \
//        } \
//        ptype##s_accumulate(sum, scratch, j); \
//        npoints -= j; \
//    } \
//}

#ifndef SCRATCH_LIMIT
# ifdef __wasm__
#  define SCRATCH_LIMIT (45 * 1024)
# else
	 /* Performance with 144K scratch is within 1-2-3% from optimal */
#  define SCRATCH_LIMIT (144 * 1024)
# endif
#endif

    //ADDITION_BTREE(blst_p1, G1Point, 384, fp, limb, get_BLST_Rx().p2)
    //    ADDITION_BTREE(blst_p2, G2Point, 384x, fp2, Fp, get_BLST_Rx().p2)



    /*
    * This implementation uses explicit addition formula:
    *
    * λ = (Y₂-Y₁)/(X₂-X₁)
    * X₃ = λ²-(X₁+X₂)
    * Y₃ = λ⋅(X₁-X₃)-Y₁
    *
    * But since we don't know if we'll have to add point to itself, we need
    * to eventually resort to corresponding doubling formula:
    *
    * λ = 3X₁²/2Y₁
    * X₃ = λ²-2X₁
    * Y₃ = λ⋅(X₁-X₃)-Y₁
    *
    * The formulae use prohibitively expensive inversion, but whenever we
    * have a lot of affine points to accumulate, we can amortize the cost
    * by applying Montgomery's batch inversion approach. As a result,
    * asymptotic[!] per-point cost for addition is as small as 5M+1S. For
    * comparison, ptype##_dadd takes 8M+5S. In practice, all things
    * considered, the improvement coefficient varies from 60% to 85%
    * depending on platform and curve.
    *
    * THIS IMPLEMENTATION IS *NOT* CONSTANT-TIME. [But if there is an
    * application that requires constant time-ness, speak up!]
    */

    /*
    * Calculate λ's numerator and denominator.
    *
    * input:	A	x1	y1	-
    *		B	x2	y2	-
    * output:
    * if A!=B:	A	x1	y1	(x2-x1)*mul_acc
    *		B	x2+x1	y2-y1	(x2-x1)
    *
    * if A==B:	A	x	y	2y*mul_acc
    *		B	2x	3*x^2	2y
    *
    * if A==-B:	A	0	0	1*mul_acc
    *		B	0	3*x^2	0
    */
    template<typename Field>
    void point_head(Point<Field> AB[2], const Field* mul_acc)
    { 
        Point<Field>& A = AB[0];
        Point<Field>& B = AB[1];
        limb_t inf = vec_is_zero(&A, sizeof(AffinePoint<Field>)) | 
                     vec_is_zero(&B, sizeof(AffinePoint<Field>));
        const Field zero = Field::zero();
    
        Field::sub(B.Z, B.X, A.X);		/* X2-X1  */ 
        Field::add(B.X, B.X, A.X);		/* X2+X1  */ 
        Field::add(A.Z, B.Y, A.Y);		/* Y2+Y1  */ 
        Field::sub(B.Y, B.Y, A.Y);		/* Y2-Y1  */ 
        if (vec_is_zero(B.Z, sizeof(B.Z))) {	/* X2==X1 */ 
            inf = vec_is_zero(&A.Z, sizeof(A.Z));	
            vec_select(&B.X, &A.Z, &B.X, sizeof(B.X), inf); 
            Field::sqr(B.Y, A.X);		
            Field::mul_by_3(B.Y, B.Y);		/* 3*X1^2 */ 
            vec_copy(&B.Z, &A.Z, sizeof(B.Z));	/* 2*Y1   */ 
        }						/* B.Y is numenator    */ 
						    /* B.Z is denominator  */ 
        vec_select(&A.X, &B.X, &A.X, sizeof(A.X), inf); 
        vec_select(&A.Y, &A.Z, &A.Y, sizeof(A.Y), inf); 
        vec_select(&A.Z, &Field::one(), &B.Z, sizeof(A.Z), inf);
        vec_select(&B.Z, &zero,& B.Z, sizeof(B.Z), inf); 
        if (mul_acc != NULL) 
            Field::mul(A.Z, A.Z, *mul_acc);	/* chain multiplication */
    }

    /*
    * Calculate λ and resulting coordinates.
    *
    * input:	A		x1			y1		-
    *		B		x2+x1			nominator	-
    * 		lambda		1/denominator
    * output:	D		x3=(nom/den)^2-(x2+x1)	y3=(nom/den)(x1-x3)-y1
    */
    template<typename Field>
    void point_tail(Point<Field>& D, Point<Field> AB[2], Field& lambda)
    { 
        Point<Field>& A = AB[0];
        Point<Field>& B = AB[1];
        Field llambda; 
        limb_t inf = vec_is_zero(&B.Z, sizeof(B.Z)); 
    
        Field::mul(lambda, lambda, B.Y);		/* λ = (Y2-Y1)/(X2-X1)  */ 
						    /* alt. 3*X1^2/2*Y1     */ 
        Field::sqr(llambda, lambda); 
        Field::sub(D.X, llambda, B.X);		/* X3 = λ^2-X1-X2       */ 
    
        Field::sub(D.Y, A.X, D.X);   
        Field::mul(D.Y, D.Y, lambda); 
        Field::sub(D.Y, D.Y, A.Y);		/* Y3 = λ*(X1-X3)-Y1    */ 
    
        vec_select(&D.X, &A.X, &D.X, 2*sizeof(D.X), inf); 
        vec_select(&B.Z, &Field::one(), &B.Z, sizeof(B.Z), inf);
    }

	/* sum = sum + Σ points[i] for i in 0..n-1.
    * 
	* This routine uses Montgomery's batch inversion to
	* speed up accumulation of multiple affine points.
    * 
    * |points[]| is volatile buffer with |X|s and |Y|s initially holding
    * input affine coordinates, and with |Z|s being used as additional
    * temporary storage [unrelated to Jacobian coordinates]. |sum| is
    * in-/output, initialize to infinity accordingly.
    */
	template<typename Field>
    void points_accumulate(Point<Field>& sum, Point<Field> points[], size_t n)
    { 
        Point<Field>* dst;
        void *mul_acc; 
        size_t i; 
    
        while (n >= 16) { 
            if (n & 1) 
                Point<Field>::dadd(sum, sum, *(const AffinePoint<Field> *)points++);
            n /= 2; 
            for (mul_acc = NULL, i = n; i--; mul_acc = points->Z, points += 2) 
                point_head<Field>(points, (const Field*)mul_acc);
    
            Field::reciprocal(points[-2].Z, points[-2].Z); /* 1/∏ Zi */
    
            for (dst = points, i = n; --i;) { 
                dst--; points -= 2; 
                Field::mul(points[-2].Z, points[0].Z, points[-2].Z); 
                point_tail<Field>(*dst, points, points[-2].Z);
                Field::mul(points[-2].Z, points[0].Z, points[1].Z); 
            } 
            dst--; points -= 2; 
            point_tail<Field>(*dst, points, points[0].Z);
            points = dst; 
        } 
        while (n--) 
            Point<Field>::dadd(sum, sum, *(const AffinePoint<Field>*)points++);
    } 
    
	template<typename Field>
    void points_add(Point<Field>& sum, const AffinePoint<Field>*const points[],
                                   size_t npoints) 
    { 
        const size_t stride = SCRATCH_LIMIT / sizeof(Point<Field>);
        Point<Field>*scratch = (Point<Field>*)malloc((npoints > stride ? stride : npoints) *
                                sizeof(Point<Field>));
        const AffinePoint<Field>*point = NULL;
    
        vec_zero(sum, sizeof(*sum)); 
        while (npoints) { 
            size_t i, j = npoints > stride ? stride : npoints; 
            for (i=0; i<j; i++) { 
                point = *points ? *points++ : point+1; 
                vec_copy(&scratch[i], point, sizeof(*point)); 
            } 
            points_accumulate<Field>(sum, scratch, j); 
            npoints -= j; 
        } 
		free(scratch);
    }
}