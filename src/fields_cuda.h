///*
// * Copyright Supranational LLC
// * Licensed under the Apache License, Version 2.0, see LICENSE for details.
// * SPDX-License-Identifier: Apache-2.0
// */
//#ifndef BLST_FIELDS_CUDA_H
//#define BLST_FIELDS_CUDA_H
//
//#include "fields.h"
//BLST_EXTERN_C_BEGIN
//
//
//// Vectorized addition of Fp elements
//void add_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num);
//
//// Vectorized subtraction of Fp elements
//void sub_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num);
//// Vectorized multiplication of Fp elements
//void mul_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num);
//
//// Vectorized squaring of Fp elements
//void sqr_fp_vec(vec384* out, const vec384* a, size_t num);
//
//// Vectorized negation of Fp elements
//void cneg_fp_vec(vec384* out, const vec384* a, int flag, size_t num);
//
//// Vectorized multiplication of Fp elements by 3
//void mul3_fp_vec(vec384* out, const vec384* a, size_t num);
//
//// Vectorized multiplication of Fp elements by 8
//BLST_INLINE void mul8_fp_vec(vec384* out, const vec384* a, size_t num);
//
//// Vectorized division of Fp elements by 2
//BLST_INLINE void div2_fp_vec(vec384* out, const vec384* a
//	, size_t num);
//
//
//// Vectorized left shift of Fp elements
//BLST_INLINE void lshift_fp_vec(vec384* out, const vec384* a
//	, size_t count, size_t num);
//
//
//// Vectorized conversion from a regular representation to a montgomery Fp
//BLST_INLINE void to_fp_vec(vec384* out, const vec384* a
//	, size_t num);
//
//
//// Vectorized conversion from a montgomery Fp to a regular representation
//BLST_INLINE void from_fp_vec(vec384* out, const vec384* a,
//	size_t num);
//
//
//BLST_EXTERN_C_END
//#endif /* BLST_FIELDS_CUDA_H */
