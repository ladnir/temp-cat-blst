/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "point.h"
#include "fields.h"
BLST_EXTERN_C_BEGIN
 /*
  * y^2 = x^3 + A'*x + B', isogenous one
  */
extern const vec384x Aprime_E2;
extern const vec384x Bprime_E2;

void map_fp2_times_Zz(vec384x map[], const vec384x isogeny_map[],
    const vec384x Zz_powers[], size_t n);
void map_fp2(vec384x acc, const vec384x x, const vec384x map[], size_t n);
void isogeny_map_to_E2(POINTonE2* out, const POINTonE2* p);
void map_to_isogenous_E2(POINTonE2* p, const vec384x u);

/*
 * As per suggestions in "7. Clearing the cofactor" at
 * https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06
 */
void POINTonE2_add_n_dbl(POINTonE2* out, const POINTonE2* p, size_t n);

void POINTonE2_times_minus_z(POINTonE2* out, const POINTonE2* in);

//void psi(POINTonE2* out, const POINTonE2* in);

void clear_cofactor(POINTonE2* out, const POINTonE2* p);
/*
 * |u|, |v| are expected to be in Montgomery representation
 */
void map_to_g2(POINTonE2* out, const vec384x u, const vec384x v);

void blst_map_to_g2(POINTonE2* out, const vec384x u, const vec384x v);

void Encode_to_G2(POINTonE2* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

void blst_encode_to_g2(POINTonE2* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

void Hash_to_G2(POINTonE2* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

void blst_hash_to_g2(POINTonE2* p, const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

bool_t POINTonE2_in_G2(const POINTonE2* P);

int blst_p2_in_g2(const POINTonE2* p);

int blst_p2_affine_in_g2(const POINTonE2_affine* p);
BLST_EXTERN_C_END
