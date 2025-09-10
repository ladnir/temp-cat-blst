#include "fields.h"

const vec384fp12* blst_fp12_one(void);
void blst_fp12_mul(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
int blst_fp12_is_one(const vec384fp12 a);
void blst_fp12_cyclotomic_sqr(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_sqr(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_conjugate(vec384fp12 a);
void blst_fp12_inverse(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_frobenius_map(vec384fp12 ret, const vec384fp12 a, size_t n);
int blst_fp12_is_equal(const vec384fp12 a, const vec384fp12 b);


void blst_fp12_from_bendian(vec384fp12 ret, const unsigned char serialized[12 * 48]);
void blst_bendian_from_fp12(unsigned char ret[48 * 12], const vec384fp12 a);
size_t blst_fp12_sizeof(void);
int blst_fp12_is_one(const vec384fp12 a);
void blst_fp12_one_copy(vec384fp12 ret);

void blst_fp12_pow(vec384fp12 result, const vec384fp12 base, const limb_t* exponent, size_t bits_in_exponent);




