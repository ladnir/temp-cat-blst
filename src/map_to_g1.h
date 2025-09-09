/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "point.h"
#include "fields.h"

 /*
  * y^2 = x^3 + A'*x + B', isogenous one
  */
extern const vec384 Aprime_E1;

extern const vec384 Bprime_E1;

void map_fp_times_Zz(vec384 map[], const vec384 isogeny_map[],
    const vec384 Zz_powers[], size_t n);

void map_fp(vec384 acc, const vec384 x, const vec384 map[], size_t n);

void isogeny_map_to_E1(POINTonE1* out, const POINTonE1* p);

void map_to_isogenous_E1(POINTonE1* p, const vec384 u);

void POINTonE1_add_n_dbl(POINTonE1* out, const POINTonE1* p, size_t n);

void POINTonE1_times_minus_z(POINTonE1* out, const POINTonE1* in);
/*
 * |u|, |v| are expected to be in Montgomery representation
 */
void map_to_g1(POINTonE1* out, const vec384 u, const vec384 v);

void blst_map_to_g1(POINTonE1* out, const vec384 u, const vec384 v);

void Encode_to_G1(POINTonE1* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);


void blst_encode_to_g1(POINTonE1* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

void Hash_to_G1(POINTonE1* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);


void blst_hash_to_g1(POINTonE1* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

void sigma(POINTonE1* out, const POINTonE1* in);

bool_t POINTonE1_in_G1(const POINTonE1* P);

int blst_p1_in_g1(const POINTonE1* p);

int blst_p1_affine_in_g1(const POINTonE1_affine* p);
