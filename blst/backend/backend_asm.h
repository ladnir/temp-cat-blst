#include "blst/config.h"
#include "blst/common/defines.h"

#ifdef BLST_ASM
BLST_EXTERN_C_BEGIN

///*
// * Assembly subroutines...
// */
//#if defined(__ADX__) /* e.g. -march=broadwell */ && !defined(__BLST_PORTABLE__)\
//                                                 && !defined(__BLST_NO_ASM__)
//# define mul_mont_sparse_256 mulx_mont_sparse_256
//# define sqr_mont_sparse_256 sqrx_mont_sparse_256
//# define from_mont_256 fromx_mont_256
//# define redc_mont_256 redcx_mont_256
//# define mul_mont_384 mulx_mont_384
//# define sqr_mont_384 sqrx_mont_384
//# define sqr_n_mul_mont_384 sqrx_n_mul_mont_384
//# define sqr_n_mul_mont_383 sqrx_n_mul_mont_383
//# define mul_384 mulx_384
//# define sqr_384 sqrx_384
//# define redc_mont_384 redcx_mont_384
//# define from_mont_384 fromx_mont_384
//# define sgn0_pty_mont_384 sgn0x_pty_mont_384
//# define sgn0_pty_mont_384x sgn0x_pty_mont_384x
//# define ct_inverse_mod_384 ctx_inverse_mod_384
//#endif
//
//#if defined(__ADX__) /* e.g. -march=broadwell */ && !defined(__BLST_PORTABLE__)
//# define mul_mont_384x mulx_mont_384x
//# define sqr_mont_384x sqrx_mont_384x
//# define sqr_mont_382x sqrx_mont_382x
//# define mul_382x mulx_382x
//# define sqr_382x sqrx_382x
//#endif


/* ret = bit ? a : b */
void vec_select_32(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_48(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_96(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_144(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_192(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_288(void* ret, const void* a, const void* b, bool_t sel_a);	// asm


// asm: is `a` a vector of zeros. requires that a is a multiple of 16 bytes in length.
bool_t vec_is_zero_16x(const void* a, size_t num);

// asm: is `a` equal to `b`. requires that each is a multiple of 16 bytes in length.
bool_t vec_is_equal_16x(const void* a, const void* b, size_t num);

// asm: 256 bit montgomery multiplication, sparse modulus version
void mulx_mont_sparse_256(vec256 ret, const vec256 a, const vec256 b,
	const vec256 p, limb_t n0);

// asm: 256 bit montgomery squaring, sparse modulus version
void sqrx_mont_sparse_256(vec256 ret, const vec256 a, const vec256 p, limb_t n0);

// asm: 256 bit montgomery squaring, sparse modulus version
//void cc_sqr_mont_sparse_256(vec256 ret, const vec256 a, const vec256 p, limb_t n0);

/* @brief asm: Performs Montgomery reduction on a 256 bit integer.
*
* Computes `ret = a / R mod p`, where `R = 2^256`. n0 is the Montgomery constant.
*/
void redcx_mont_256(vec256 ret, const vec512 a, const vec256 p, limb_t n0);

/**
 * @brief asm: Converts a large integer from Montgomery form.
 *
 * Computes `a * R^(-1) mod p` to convert `a` from Montgomery form to standard form.
 * The result is stored in `ret`. n0 is the Montgomery constant.
 */
void fromx_mont_256(vec256 ret, const vec256 a, const vec256 p, limb_t n0);


/* @brief asm: Adds two 256-bit integers modulo a given modulus.
*/
void add_mod_256(vec256 ret, const vec256 a, const vec256 b, const vec256 p);


/* @brief asm: Subtracts two 256-bit integers modulo a given modulus.
*/
void sub_mod_256(vec256 ret, const vec256 a, const vec256 b, const vec256 p);

/* @brief asm: Multiplies a 256-bit integer by 3 modulo a given modulus.
*/
void mul_by_3_mod_256(vec256 ret, const vec256 a, const vec256 p);

/* @brief asm: Conditionally negates a 256-bit integer modulo a given modulus.
*/
void cneg_mod_256(vec256 ret, const vec256 a, bool_t flag, const vec256 p);

/* @brief asm: Left shifts a 256-bit integer modulo a given modulus.
*/
void lshift_mod_256(vec256 ret, const vec256 a, size_t count, const vec256 p);

/* @brief asm: Right shifts a 256-bit integer modulo a given modulus.
*/
void rshift_mod_256(vec256 ret, const vec256 a, size_t count, const vec256 p);

/* @brief asm: Checks if a 256-bit integer is less than a given modulus.
*/
limb_t check_mod_256(const pow256 a, const vec256 p);

/* @brief asm: Adds two 256-bit integers modulo a given modulus and checks if the result is zero.
*/
limb_t add_n_check_mod_256(pow256 ret, const pow256 a, const pow256 b, const vec256 p);

/* @brief asm: Subtracts two 256-bit integers modulo a given modulus and checks if the result is zero.
*/
limb_t sub_n_check_mod_256(pow256 ret, const pow256 a, const pow256 b, const vec256 p);


/* @brief asm: 384 bit montgomery multiplication mod p. n0 is the montgomery constant. */ 
void mulx_mont_384(vec384 ret, const vec384 a, const vec384 b, const vec384 p, limb_t n0);

/* @brief asm: 384 bit montgomery squaring mod p. n0 is the montgomery constant. */
void sqrx_mont_384(vec384 ret, const vec384 a, const vec384 p, limb_t n0);

///* @brief asm: Squares a 384-bit integer multiple times and multiplies the result in Montgomery form. */
//void sqrx_n_mul_mont_384(vec384 ret, const vec384 a, size_t count,
//	const vec384 p, limb_t n0, const vec384 b);

/* @brief asm: 382 bit montgomery repeated squaring (a^(2^count) * b * R^(-1)) mod p. n0 is the montgomery constant. 
 *
*"383" - bit variant omits full reductions at the ends of squarings,
* which results in up to ~15 % improvement.[One can improve further
* by omitting full reductions even after multiplications and
*performing final reduction at the very end of the chain.]
* 
* see sqr_n_mul_fp
*/
void sqrx_n_mul_mont_383(vec384 ret, const vec384 a, size_t count,
	const vec384 p, limb_t n0, const vec384 b);

#if defined(BLST_ADX) || defined(__MSC_VER)
/* @brief asm: Performs Montgomery reduction on a 384 bit integer. n0 is the montgomery constant. */
void redcx_mont_384(vec384 ret, const vec768 a, const vec384 p, limb_t n0);
#endif

#if !defined(BLST_ADX) || defined(__MSC_VER)
/* @brief asm: Performs Montgomery reduction on a 384 bit integer. n0 is the montgomery constant.  */
void redc_mont_384(vec384 ret, const vec768 a, const vec384 p, limb_t n0);
#endif

/* @brief asm: Converts a large integer from Montgomery form. n0 is the montgomery constant. */
void fromx_mont_384(vec384 ret, const vec384 a, const vec384 p, limb_t n0);

/* @brief asm: Computes the sign and parity of a 384-bit integer in Montgomery form modulo a given modulus. */
limb_t sgn0x_pty_mont_384(const vec384 a, const vec384 p, limb_t n0);

/* @brief asm: Computes the sign and parity of a 384-bit complex number in Montgomery form modulo a given modulus. */
limb_t sgn0x_pty_mont_384x(const vec384x a, const vec384 p, limb_t n0);

/* @brief asm: Compute sign and parity of a in Montgomery form mod p. */
limb_t sgn0_pty_mod_384(const vec384 a, const vec384 p);

/* @brief asm: Compute sign and parity of a complex number mod p. */
limb_t sgn0_pty_mod_384x(const vec384x a, const vec384 p);

/* @brief asm: Conditionally negates a 384-bit integer modulo a given modulus. */
void add_mod_384(vec384 ret, const vec384 a, const vec384 b, const vec384 p);

/* @brief asm: Conditionally negates a 384-bit integer modulo a given modulus. */
void sub_mod_384(vec384 ret, const vec384 a, const vec384 b, const vec384 p);

/* @brief asm: Multiplies a 384-bit integer by 8 modulo a given modulus. */
void mul_by_8_mod_384(vec384 ret, const vec384 a, const vec384 p);

/* @brief asm: Multiplies a 384-bit integer by 3 modulo a given modulus. */
void mul_by_3_mod_384(vec384 ret, const vec384 a, const vec384 p);

/* @brief asm: Conditionally negates a 384-bit integer modulo a given modulus. */
void cneg_mod_384(vec384 ret, const vec384 a, bool_t flag, const vec384 p);

/* @brief asm: right shifts a 384-bit integer modulo a given modulus. */
void rshift_mod_384(vec384 ret, const vec384 a, size_t count, const vec384 p);

/* @brief asm: divides a 384-bit integer by 2 modulo a given modulus. */
void div_by_2_mod_384(vec384 ret, const vec384 a, const vec384 p);

/* @brief asm: Computes the modular inverse of a 384-bit integer modulo a given modulus. 
* modx should be `-p^(-1) mod R` where `R=2^(n*BLST_LIMB_SIZE)`.
*/
void ctx_inverse_mod_384(vec768 ret, const vec384 inp, const vec384 mod, const vec384 modx);

/* @brief asm: Computes the modular inverse of a 256-bit integer modulo a given modulus.
* modx should be `-p^(-1) mod R` where `R=2^(n*BLST_LIMB_SIZE)`.
*/
void ct_inverse_mod_256(vec512 ret, const vec256 inp, const vec256 mod, const vec256 modx);

/* @brief asm: Checks if a 384-bit integer is a quadratic residue modulo a given modulus. */
bool_t ct_is_square_mod_384(const vec384 inp, const vec384 mod);

/* @brief asm: Performs Montgomery multiplication on two 384-bit complex numbers. */
void mulx_mont_384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p, limb_t n0);

/* @brief asm: Squares a 384-bit complex number in Montgomery form . */
void sqrx_mont_384x(vec384x ret, const vec384x a, const vec384 p, limb_t n0);

/* @brief asm: Squares a 382-bit complex in Montgomery form. */
void sqrx_mont_382x(vec384x ret, const vec384x a, const vec384 p, limb_t n0);

/* @brief asm: adds two 384-bit complex integers modulo a given modulus. */
void add_mod_384x(vec384x ret, const vec384x a, const vec384x b, const vec384 p);

/* @brief asm: subtracts two 384-bit complex integers modulo a given modulus. */
void sub_mod_384x(vec384x ret, const vec384x a, const vec384x b, const vec384 p);

/* @brief asm: multiple complex by 8 mod p */
void mul_by_8_mod_384x(vec384x ret, const vec384x a, const vec384 p);

/* @brief asm: multiple complex by 3 mod p */
void mul_by_3_mod_384x(vec384x ret, const vec384x a, const vec384 p);

/* @brief asm: multiply complex by i */
void mul_by_1_plus_i_mod_384x(vec384x ret, const vec384x a, const vec384 p);

/* @brief asm: sha256 block processing */
void blst_sha256_block_data_order(unsigned int* v, const void* inp, size_t blocks);

/* @brief asm: copy sha256 hash state */
void blst_sha256_hcopy(unsigned int dst[8], const unsigned int src[8]);

/* @brief asm: emit sha256 hash */
void blst_sha256_emit(unsigned char md[32], const unsigned int h[8]);

/* @brief asm: copy memory */
void blst_sha256_bcopy(void* dst_, const void* src_, size_t len);

/* @brief asm: Divides a 3 - limb dividend by a 2 - limb divisor.*/
limb_t div_3_limbs(const limb_t dividend_top[2], limb_t d_lo, limb_t d_hi);

/* @brief asm: Computes the quotient and remainder of 128 bit integer division.*/
limb_t quot_rem_128(limb_t* quot_rem, const limb_t* divisor, limb_t quotient);

/* @brief Computes the quotient and remainder of a 64 bit integer division. */
limb_t quot_rem_64(limb_t* quot_rem, const limb_t* divisor, limb_t quotient);



void add_mod_384x384(vec768 ret, const vec768 a, const vec768 b,
	const vec384 p);

void sub_mod_384x384(vec768 ret, const vec768 a, const vec768 b,
	const vec384 p);


#ifdef BLST_ADX

void mulx_382x(vec768x ret, const vec384x a, const vec384x b, const vec384 p);   /* +~6% in Miller loop */

void sqrx_382x(vec768x ret, const vec384x a, const vec384 p);      /* +~5% in final exponentiation */

#else

void mul_382x(vec768x ret, const vec384x a, const vec384x b, const vec384 p);   /* +~6% in Miller loop */

void sqr_382x(vec768x ret, const vec384x a, const vec384 p);      /* +~5% in final exponentiation */

#endif // BLST_ADX

BLST_EXTERN_C_END
#endif // BLST_ASM