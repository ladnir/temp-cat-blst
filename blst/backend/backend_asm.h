#include "blst/config.h"
#include "blst/common/defines.h"


#ifdef BLST_ASM
BLST_EXTERN_C_BEGIN_


/* ret = bit ? a : b */
void vec_select_32(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
void vec_select_48(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
void vec_select_96(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
void vec_select_144(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
void vec_select_192(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
void vec_select_288(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm


#define asm_pow256  /*pow256*/ void*
#define asm_Fr      /*Fr    */ void*
#define asm_UFr     /*UFr   */ void*
#define asm_Fp      /*Fp    */ void*
#define asm_UFp     /*UFp   */ void*
#define asm_Fp2     /*Fp2   */ void*
#define asm_UFp2    /*UFp2  */ void*

/* @brief asm: Adds two 256-bit integers modulo a given modulus and checks if the result is zero.
*/
blst::limb_t add_n_check_mod_256(asm_pow256 ret, const asm_pow256 a, const asm_pow256 b, const asm_Fr p);

/* @brief asm: Subtracts two 256-bit integers modulo a given modulus and checks if the result is zero.
*/
blst::limb_t sub_n_check_mod_256(asm_pow256 ret, const asm_pow256 a, const asm_pow256 b, const asm_Fr p);


// asm: 256 bit montgomery multiplication, sparse modulus version
void mulx_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr b,
    const asm_Fr p, blst::limb_t n0);

// asm: 256 bit montgomery squaring, sparse modulus version
void sqrx_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr p, blst::limb_t n0);

// asm: 256 bit montgomery squaring, sparse modulus version
//void cc_sqr_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr p, blst::limb_t n0);

/* @brief asm: Performs Montgomery reduction on a 256 bit integer.
*
* Computes `ret = a / R mod p`, where `R = 2^256`. n0 is the Montgomery constant.
*/
void redcx_mont_256(asm_Fr ret, const asm_UFr a, const asm_Fr p, blst::limb_t n0);

/**
 * @brief asm: Converts a large integer from Montgomery form.
 *
 * Computes `a * R^(-1) mod p` to convert `a` from Montgomery form to standard form.
 * The result is stored in `ret`. n0 is the Montgomery constant.
 */
void fromx_mont_256(asm_Fr ret, const asm_Fr a, const asm_Fr p, blst::limb_t n0);


/* @brief asm: Adds two 256-bit integers modulo a given modulus.
*/
void add_mod_256(asm_Fr ret, const asm_Fr a, const asm_Fr b, const asm_Fr p);


/* @brief asm: Subtracts two 256-bit integers modulo a given modulus.
*/
void sub_mod_256(asm_Fr ret, const asm_Fr a, const asm_Fr b, const asm_Fr p);

/* @brief asm: Multiplies a 256-bit integer by 3 modulo a given modulus.
*/
void mul_by_3_mod_256(asm_Fr ret, const asm_Fr a, const asm_Fr p);

/* @brief asm: Conditionally negates a 256-bit integer modulo a given modulus.
*/
void cneg_mod_256(asm_Fr ret, const asm_Fr a, blst::bool_t flag, const asm_Fr p);

/* @brief asm: Left shifts a 256-bit integer modulo a given modulus.
*/
void lshift_mod_256(asm_Fr ret, const asm_Fr a, size_t count, const asm_Fr p);

/* @brief asm: Right shifts a 256-bit integer modulo a given modulus.
*/
void rshift_mod_256(asm_Fr ret, const asm_Fr a, size_t count, const asm_Fr p);

/* @brief asm: Checks if a 256-bit integer is less than a given modulus.
*/
blst::limb_t check_mod_256(const asm_pow256 a, const asm_Fr p);

// asm: is `a` a vector of zeros. requires that a is a multiple of 16 bytes in length.
blst::bool_t vec_is_zero_16x(const void* a, size_t num);

// asm: is `a` equal to `b`. requires that each is a multiple of 16 bytes in length.
blst::bool_t vec_is_equal_16x(const void* a, const void* b, size_t num);


/* @brief asm: 384 bit montgomery multiplication mod p. n0 is the montgomery constant. */ 
void mulx_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp b, const asm_Fp p, blst::limb_t n0);

/* @brief asm: 384 bit montgomery squaring mod p. n0 is the montgomery constant. */
void sqrx_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p, blst::limb_t n0);

///* @brief asm: Squares a 384-bit integer multiple times and multiplies the result in Montgomery form. */
//void sqrx_n_mul_mont_384(asm_Fp ret, const asm_Fp a, size_t count,
//	const asm_Fp p, blst::limb_t n0, const asm_Fp b);

/* @brief asm: 382 bit montgomery repeated squaring (a^(2^count) * b * R^(-1)) mod p. n0 is the montgomery constant. 
 *
*"383" - bit variant omits full reductions at the ends of squarings,
* which results in up to ~15 % improvement.[One can improve further
* by omitting full reductions even after multiplications and
*performing final reduction at the very end of the chain.]
* 
* see sqr_n_Fp::mul
*/
void sqrx_n_mul_mont_383(asm_Fp ret, const asm_Fp a, size_t count,
	const asm_Fp p, blst::limb_t n0, const asm_Fp b);

#if defined(BLST_ADX) || defined(__MSC_VER)
/* @brief asm: Performs Montgomery reduction on a 384 bit integer. n0 is the montgomery constant. */
void redcx_mont_384(asm_Fp ret, const asm_UFp a, const asm_Fp p, blst::limb_t n0);
#endif

#if !defined(BLST_ADX) || defined(__MSC_VER)
/* @brief asm: Performs Montgomery reduction on a 384 bit integer. n0 is the montgomery constant.  */
void redc_mont_384(asm_Fp ret, const asm_UFp a, const asm_Fp p, blst::limb_t n0);
#endif

/* @brief asm: Converts a large integer from Montgomery form. n0 is the montgomery constant. */
void fromx_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p, blst::limb_t n0);

/* @brief asm: Computes the sign and parity of a 384-bit integer in Montgomery form modulo a given modulus. */
blst::limb_t sgn0x_pty_mont_384(const asm_Fp a, const asm_Fp p, blst::limb_t n0);

/* @brief asm: Computes the sign and parity of a 384-bit complex number in Montgomery form modulo a given modulus. */
blst::limb_t sgn0x_pty_mont_384x(const asm_Fp2 a, const asm_Fp p, blst::limb_t n0);

/* @brief asm: Compute sign and parity of a in Montgomery form mod p. */
blst::limb_t sgn0_pty_mod_384(const asm_Fp a, const asm_Fp p);

/* @brief asm: Compute sign and parity of a complex number mod p. */
blst::limb_t sgn0_pty_mod_384x(const asm_Fp2 a, const asm_Fp p);

/* @brief asm: Conditionally negates a 384-bit integer modulo a given modulus. */
void add_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp b, const asm_Fp p);

/* @brief asm: Conditionally negates a 384-bit integer modulo a given modulus. */
void sub_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp b, const asm_Fp p);

/* @brief asm: Multiplies a 384-bit integer by 8 modulo a given modulus. */
void mul_by_8_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp p);

/* @brief asm: Multiplies a 384-bit integer by 3 modulo a given modulus. */
void mul_by_3_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp p);

/* @brief asm: Conditionally negates a 384-bit integer modulo a given modulus. */
void cneg_mod_384(asm_Fp ret, const asm_Fp a, blst::bool_t flag, const asm_Fp p);

/* @brief asm: right shifts a 384-bit integer modulo a given modulus. */
void rshift_mod_384(asm_Fp ret, const asm_Fp a, size_t count, const asm_Fp p);

/* @brief asm: left shifts a 384-bit integer modulo a given modulus. */
void lshift_mod_384(asm_Fp ret, const asm_Fp a, size_t count, const asm_Fp p);

/* @brief asm: divides a 384-bit integer by 2 modulo a given modulus. */
void div_by_2_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp p);

/* @brief asm: Computes the modular inverse of a 384-bit integer modulo a given modulus. 
* modx should be `-p^(-1) mod R` where `R=2^(n*BLST_LIMB_SIZE)`.
*/
void ctx_inverse_mod_384(asm_UFp ret, const asm_Fp inp, const asm_Fp mod, const asm_Fp modx);

/* @brief asm: Computes the modular inverse of a 256-bit integer modulo a given modulus.
* modx should be `-p^(-1) mod R` where `R=2^(n*BLST_LIMB_SIZE)`.
*/
void ct_inverse_mod_256(asm_UFr ret, const asm_Fr inp, const asm_Fr mod, const asm_Fr modx);

/* @brief asm: Checks if a 384-bit integer is a quadratic residue modulo a given modulus. */
blst::bool_t ct_is_square_mod_384(const asm_Fp inp, const asm_Fp mod);

/* @brief asm: Performs Montgomery multiplication on two 384-bit complex numbers. */
void mulx_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b,
	const asm_Fp p, blst::limb_t n0);

/* @brief asm: Squares a 384-bit complex number in Montgomery form . */
void sqrx_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p, blst::limb_t n0);

/* @brief asm: Squares a 382-bit complex in Montgomery form. */
void sqrx_mont_382x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p, blst::limb_t n0);

/* @brief asm: adds two 384-bit complex integers modulo a given modulus. */
void add_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b, const asm_Fp p);

/* @brief asm: subtracts two 384-bit complex integers modulo a given modulus. */
void sub_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b, const asm_Fp p);

/* @brief asm: multiple complex by 8 mod p */
void mul_by_8_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p);

/* @brief asm: multiple complex by 3 mod p */
void mul_by_3_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p);

/* @brief asm: multiply complex by i */
void mul_by_1_plus_i_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p);

#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
void mul_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr b,
    const asm_Fr p, blst::limb_t n0);
void sqr_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr p, blst::limb_t n0);
void redc_mont_256(asm_Fr ret, const asm_UFr a, const asm_Fr p, blst::limb_t n0);
void from_mont_256(asm_Fr ret, const asm_Fr a, const asm_Fr p, blst::limb_t n0);

void mul_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp b, const asm_Fp p, blst::limb_t n0);
void sqr_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p, blst::limb_t n0);
void from_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p, blst::limb_t n0);
void sqr_n_mul_mont_383(asm_Fp ret, const asm_Fp a, size_t count,
    const asm_Fp p, blst::limb_t n0, const asm_Fp b);
blst::limb_t sgn0_pty_mont_384(const asm_Fp a, const asm_Fp p, blst::limb_t n0);
void ct_inverse_mod_384(asm_UFp ret, const asm_Fp inp, const asm_Fp mod, const asm_Fp modx);

void mul_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p, blst::limb_t n0);
void sqr_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p, blst::limb_t n0);

#define mulx_mont_sparse_256 mul_mont_sparse_256
#define sqrx_mont_sparse_256 sqr_mont_sparse_256
#define redcx_mont_256 redc_mont_256
#define fromx_mont_256 from_mont_256

#define mulx_mont_384 mul_mont_384
#define sqrx_mont_384 sqr_mont_384
#define fromx_mont_384 from_mont_384
#define sqrx_n_mul_mont_383 sqr_n_mul_mont_383
#define sgn0x_pty_mont_384 sgn0_pty_mont_384
#define ctx_inverse_mod_384 ct_inverse_mod_384

#define mulx_mont_384x mul_mont_384x
#define sqrx_mont_384x sqr_mont_384x
#endif

/* @brief asm: sha256 block processing */
void blst_sha256_block_data_order(unsigned int* v, const void* inp, size_t blocks);

/* @brief asm: copy sha256 hash state */
void blst_sha256_hcopy(unsigned int dst[8], const unsigned int src[8]);

/* @brief asm: emit sha256 hash */
void blst_sha256_emit(unsigned char md[32], const unsigned int h[8]);

/* @brief asm: copy memory */
void blst_sha256_bcopy(void* dst_, const void* src_, size_t len);

/* @brief asm: Divides a 3 - limb_t dividend by a 2 - limb_t divisor.*/
blst::limb_t div_3_limbs(const blst::limb_t dividend_top[2], blst::limb_t d_lo, blst::limb_t d_hi);

/* @brief asm: Computes the quotient and remainder of 128 bit integer division.*/
blst::limb_t quot_rem_128(blst::limb_t* quot_rem, const blst::limb_t* divisor, blst::limb_t quotient);

/* @brief Computes the quotient and remainder of a 64 bit integer division. */
blst::limb_t quot_rem_64(blst::limb_t* quot_rem, const blst::limb_t* divisor, blst::limb_t quotient);



void add_mod_384x384(asm_UFp ret, const asm_UFp a, const asm_UFp b,
	const asm_Fp p);

void sub_mod_384x384(asm_UFp ret, const asm_UFp a, const asm_UFp b,
	const asm_Fp p);


#ifdef BLST_ADX

void mulx_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp2 b, const asm_Fp p);   /* +~6% in Miller loop */

void sqrx_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp p);      /* +~5% in final exponentiation */

#else

void mul_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp2 b, const asm_Fp p);   /* +~6% in Miller loop */

void sqr_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp p);      /* +~5% in final exponentiation */

#endif // BLST_ADX


#if (defined(__x86_64__) || defined(__x86_64) || defined(_M_X64)) && \
     defined(__SHA__) /* -msha */ 
void blst_sha256_block_data_order_shaext(unsigned int* v, const void* inp, size_t blocks);
#endif



#undef asm_pow256  
#undef asm_Fr      
#undef asm_UFr     
#undef asm_Fp      
#undef asm_UFp     
#undef asm_Fp2     
#undef asm_UFp2    

BLST_EXTERN_C_END_
#endif // BLST_ASM
