/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_POINT_H__
#define __BLS12_381_ASM_POINT_H__


#include "config.h"
#include "vect.h"
#include "bytes.h"
#include "consts.h"


BLST_EXTERN_C_BEGIN

#define DECLARE_POINT(ptype, bits)                                           \
typedef struct { vec##bits X,Y,Z; } ptype;                                   \
typedef struct { vec##bits X,Y; } ptype##_affine;                            \
                                                                             \
HOST_DEVICE_FUNC void ptype##_dadd(ptype *out, const ptype *p1, const ptype *p2,	             \
                         const vec##bits a4);				                 \
HOST_DEVICE_FUNC void ptype##_dadd_affine(ptype *out, const ptype *p1,		                 \
                                            const ptype##_affine *p2);	     \
HOST_DEVICE_FUNC void ptype##_add(ptype *out, const ptype *p1, const ptype *p2);	             \
HOST_DEVICE_FUNC void ptype##_add_affine(ptype *out, const ptype *p1,		                 \
                                           const ptype##_affine *p2);	     \
HOST_DEVICE_FUNC void ptype##_double(ptype *out, const ptype *p1);		                     \
HOST_DEVICE_FUNC void ptype##_mult_w5(ptype *out, const ptype *point,		                 \
                            const byte *scalar, size_t nbits);		         \
HOST_DEVICE_FUNC void ptype##_cneg(ptype *p, limb_t cbit);			                         \
HOST_DEVICE_FUNC void ptype##_to_affine(ptype##_affine *out, const ptype *in);	             \
HOST_DEVICE_FUNC void ptype##_from_Jacobian(ptype *out, const ptype *in);		             \
                                                                             \
HOST_DEVICE_FUNC BLST_INLINE void ptype##_cswap(ptype *restrict a,			                     \
                                 ptype *restrict b, bool_t cbit) {	         \
    vec_cswap(a, b, sizeof(ptype), cbit);				                     \
}                                                                            \
HOST_DEVICE_FUNC BLST_INLINE void ptype##_ccopy(ptype *restrict a,			                     \
                                 const ptype *restrict b, bool_t cbit) {     \
    vec_select(a, b, a, sizeof(ptype), cbit);				                 \
}

#define DECLARE_PRIVATE_POINTXZ(ptype, bits)                                 \
typedef struct { vec##bits X,Z; } ptype##xz;                                 \
                                                                             \
HOST_DEVICE_FUNC void ptype##xz_ladder_pre(ptype##xz *out, const ptype *in);	                 \
HOST_DEVICE_FUNC void ptype##xz_ladder_step(ptype##xz *r, ptype##xz *s,		                 \
                                  const ptype##xz *p);			             \
HOST_DEVICE_FUNC void ptype##xz_ladder_post(ptype *ret,				                         \
                                  const ptype##xz *r, const ptype##xz *s,    \
                                  const ptype##xz *p, const vec##bits Y1);   \
                                                                             \
HOST_DEVICE_FUNC BLST_INLINE void ptype##xz_cswap(ptype##xz *restrict a,		                     \
                                   ptype##xz *restrict b, bool_t cbit) {     \
    vec_cswap(a, b, sizeof(ptype##xz), cbit);				                 \
}

DECLARE_POINT(POINTonE1, 384)

DECLARE_POINT(POINTonE2, 384x)

#ifdef __GNUC__
# pragma GCC diagnostic ignored "-Wunused-function"
#endif



#define BLST_B_E1_LIMBS \
    TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),\
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),\
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e)

#define BLST_G1_LIMBS                                                     \
 /* generator point [in Montgomery] */                                    \
    /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905  */            \
    /*    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */ \
{                                                                         \
    TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),         \
        TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),     \
        TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75)      \
},                                                                        \
/* (0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af6             */     \
/*    00db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1 << 384) % P */     \
{ TO_LIMB_T(0xbaac93d50ce72271), TO_LIMB_T(0x8c22631a7918fd8e),           \
  TO_LIMB_T(0xdd595f13570725ce), TO_LIMB_T(0x51ac582950405194),           \
  TO_LIMB_T(0x0e1c8c3fad0059c0), TO_LIMB_T(0x0bbc3efc5008a26a) },         \
{ ONE_MONT_P_LIMBS }


#define BLST_NEG_G1_LIMBS                                                    \
/* negative generator [in Montgomery] */                                     \
    /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905*/                 \
     /*    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */   \
{                                                                            \
    TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),            \
        TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),        \
        TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75)         \
},                                                                           \
/* (0x114d1d6855d545a8aa7d76c8cf2e21f267816aef1db507c9*/                     \
/*    6655b9d5caac42364e6f38ba0ecb751bad54dcd6b939c2ca << 384) % P */        \
{ TO_LIMB_T(0xff526c2af318883a), TO_LIMB_T(0x92899ce4383b0270),              \
  TO_LIMB_T(0x89d7738d9fa9d055), TO_LIMB_T(0x12caf35ba344c12a),              \
  TO_LIMB_T(0x3cff1b76964b5317), TO_LIMB_T(0x0e44d2ede9774430) },            \
{ ONE_MONT_P_LIMBS }




/*
 * y^2 = x^3 + B
 */
#define BLST_B_E2_LIMBS                                               \
  { TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),     \
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),     \
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e) },   \
  { TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),     \
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),     \
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e) }

 /* generator point [in Montgomery] */
#define BLST_G2_LIMBS                                                      \
{ /* (0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02            */     \
/*      b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 << 384) % P */    \
  { TO_LIMB_T(0xf5f28fa202940a10), TO_LIMB_T(0xb3f5fb2687b4961a),          \
    TO_LIMB_T(0xa1a893b53e2ae580), TO_LIMB_T(0x9894999d1a3caee9),          \
    TO_LIMB_T(0x6f67b7631863366b), TO_LIMB_T(0x058191924350bcd7) },        \
    /* (0x13e02b6052719f607dacd3a088274f65596bd0d09920b61a            */   \
    /*  b5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e << 384) % P */    \
    { TO_LIMB_T(0xa5a9c0759e23f606), TO_LIMB_T(0xaaa0c59dbccd60c3),        \
      TO_LIMB_T(0x3bb17e18e2867806), TO_LIMB_T(0x1b1ab6cc8541b367),        \
      TO_LIMB_T(0xc2b6ed0ef2158547), TO_LIMB_T(0x11922a097360edf3) }       \
  },                                                                       \
  { /* (0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a7            */   \
    /* 6d429a695160d12c923ac9cc3baca289e193548608b82801 << 384) % P */     \
    { TO_LIMB_T(0x4c730af860494c4a), TO_LIMB_T(0x597cfa1f5e369c5a),        \
      TO_LIMB_T(0xe7e6856caa0a635a), TO_LIMB_T(0xbbefb5e96e0d495f),        \
      TO_LIMB_T(0x07d3a975f0ef25a2), TO_LIMB_T(0x0083fd8e7e80dae5) },      \
      /* (0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af          */   \
    /*    267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be << 384) % P */  \
      { TO_LIMB_T(0xadc0fc92df64b05d), TO_LIMB_T(0x18aa270a2b1461dc),      \
        TO_LIMB_T(0x86adac6a3be4eba0), TO_LIMB_T(0x79495c4ec93da33a),      \
        TO_LIMB_T(0xe7175850a43ccaed), TO_LIMB_T(0x0b2bc2a163de1bf2) },    \
    },                                                                     \
    { { ONE_MONT_P_LIMBS }, { 0 } }

/* negative generator [in Montgomery] */
#define BLST_NEG_G2_LIMBS                                                   \
{ /* (0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02             */     \
  /*    b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 << 384) % P  */    \
  { TO_LIMB_T(0xf5f28fa202940a10), TO_LIMB_T(0xb3f5fb2687b4961a),           \
    TO_LIMB_T(0xa1a893b53e2ae580), TO_LIMB_T(0x9894999d1a3caee9),           \
    TO_LIMB_T(0x6f67b7631863366b), TO_LIMB_T(0x058191924350bcd7) },         \
    /* (0x13e02b6052719f607dacd3a088274f65596bd0d09920b61a             */   \
    /*    b5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e << 384) % P */   \
    { TO_LIMB_T(0xa5a9c0759e23f606), TO_LIMB_T(0xaaa0c59dbccd60c3),         \
      TO_LIMB_T(0x3bb17e18e2867806), TO_LIMB_T(0x1b1ab6cc8541b367),         \
      TO_LIMB_T(0xc2b6ed0ef2158547), TO_LIMB_T(0x11922a097360edf3) }        \
  },                                                                        \
  { /* (0x0d1b3cc2c7027888be51d9ef691d77bcb679afda66c73f17             */   \
    /*    f9ee3837a55024f78c71363275a75d75d86bab79f74782aa << 384) % P */   \
    { TO_LIMB_T(0x6d8bf5079fb65e61), TO_LIMB_T(0xc52f05df531d63a5),         \
      TO_LIMB_T(0x7f4a4d344ca692c9), TO_LIMB_T(0xa887959b8577c95f),         \
      TO_LIMB_T(0x4347fe40525c8734), TO_LIMB_T(0x197d145bbaff0bb5) },       \
      /* (0x13fa4d4a0ad8b1ce186ed5061789213d993923066dddaf10             */ \
      /*    40bc3ff59f825c78df74f2d75467e25e0f55f8a00fa030ed << 384) % P */ \
      { TO_LIMB_T(0x0c3e036d209afa4e), TO_LIMB_T(0x0601d8f4863f9e23),       \
        TO_LIMB_T(0xe0832636bacc0a84), TO_LIMB_T(0xeb2def362a476f84),       \
        TO_LIMB_T(0x64044f659f0ee1e9), TO_LIMB_T(0x0ed54f48d5a1caa7) }      \
    },                                                                      \
    { { ONE_MONT_P_LIMBS }, { 0 } }                                         


/* such that beta^3 - 1 = 0  */
#define BLST_beta_LIMBS                                                    \
    /* -1/2 * (1 + sqrt(-3)) = ((P-2)^(P-2)) * (1 + (P-3)^((P+1)/4)) */    \
    /* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4*/               \
    /*    897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac << 384) % P */  \
    TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),          \
    TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),          \
    TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)



/* y^2 = x^3 + B,   (4 << 384) % P */
static const vec384 BLST_B_E1_HOST = { BLST_B_E1_LIMBS };
static const POINTonE1 BLST_G1_HOST = { BLST_G1_LIMBS };
static const POINTonE1 BLST_NEG_G1_HOST = { BLST_NEG_G1_LIMBS };
static const vec384 BLST_beta_HOST = { BLST_beta_LIMBS };


/*
 * y^2 = x^3 + B, 4 + 4*i 
 */
const vec384x BLST_B_E2_HOST = { BLST_B_E2_LIMBS };
const POINTonE2 BLST_G2_HOST = { BLST_G2_LIMBS };
const POINTonE2 BLST_NEG_G2_HOST = { BLST_NEG_G2_LIMBS };


#if defined(BLST_CUDA) && defined(__CUDACC__)

static __constant__ vec384 BLST_B_E1_DEV = { BLST_B_E1_LIMBS };
static __constant__ POINTonE1 BLST_G1_DEV = { BLST_G1_LIMBS };
static __constant__ POINTonE1 BLST_NEG_G1_DEV = { BLST_NEG_G1_LIMBS };
static __constant__ vec384x BLST_B_E2_DEV = { BLST_B_E2_LIMBS };
static __constant__ POINTonE2 BLST_G2_DEV = { BLST_G2_LIMBS };
static __constant__ POINTonE2 BLST_NEG_G2_DEV = { BLST_NEG_G2_LIMBS };
static __constant__ vec384 BLST_beta_DEV = { BLST_beta_LIMBS };

#endif // __CUDA_ARCH__


// return the curve parameter B for the G1 curve: y^2 = x^3 + B
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec384* get_BLST_B_E1(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
    return &BLST_B_E1_DEV;
#else
    return &BLST_B_E1_HOST;
#endif
}

// return the generator point of G1 in affine coordinates, in montgomery form.
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const POINTonE1* get_BLST_G1(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
    return &BLST_G1_DEV;
#else
    return &BLST_G1_HOST;
#endif
}

// return the negative of the generator point of G1 in affine coordinates, in montgomery form.
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const POINTonE1* get_BLST_NEG_G1(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
    return &BLST_NEG_G1_DEV;
#else
    return &BLST_NEG_G1_HOST;
#endif
}

// return the curve parameter B for the G2 curve: y^2 = x^3 + B
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec384x* get_BLST_B_E2(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_B_E2_DEV;
#else
	return &BLST_B_E2_HOST;
#endif
}

// return the generator point of G2 in affine coordinates, in montgomery form.
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const POINTonE2* get_BLST_G2(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_G2_DEV;
#else
	return &BLST_G2_HOST;
#endif
}

// return the negative of the generator point of G2 in affine coordinates, in montgomery form.
HOST_DEVICE_FUNC BLST_FORCEINLINE const POINTonE2* get_BLST_NEG_G2(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_NEG_G2_DEV;
#else
	return &BLST_NEG_G2_HOST;
#endif
}


/* such that beta^3 - 1 = 0  */
/* -1/2 * (1 + sqrt(-3)) = ((P-2)^(P-2)) * (1 + (P-3)^((P+1)/4)) */
/* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
      897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac << 384) % P */
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec384* get_BLST_beta(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
    return &BLST_beta_DEV;
#else
    return &BLST_beta_HOST;
#endif
}




BLST_EXTERN_C_END
#endif
