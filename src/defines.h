#ifndef __BLS12_381_DEFINES_H__
#define __BLS12_381_DEFINES_H__

#include "config.h"

#if BLST_LIMB_SIZE == 64
typedef unsigned long long limb_t;
#elif BLST_LIMB_SIZE == 32
typedef unsigned int limb_t;
#else
# error "Unsupported BLST_LIMB_SIZE"
#endif

typedef unsigned char byte;
/*
 * Why isn't BLST_LIMB_SIZE defined as 8*sizeof(limb_t)? Because pre-processor
 * knows nothing about sizeof(anything)...
 */
#if BLST_LIMB_SIZE == 64
# define TO_LIMB_T(limb64)     limb64
#else
# define TO_LIMB_T(limb64)     (limb_t)limb64,(limb_t)(limb64>>32)
#endif

#define NLIMBS(bits)   (bits/BLST_LIMB_SIZE)

typedef limb_t vec256[NLIMBS(256)];
typedef limb_t vec512[NLIMBS(512)];
typedef limb_t vec384[NLIMBS(384)];
typedef limb_t vec768[NLIMBS(768)];
typedef vec384 vec384x[2];      /* 0 is "real" part, 1 is "imaginary" */

#define TO_BYTES(limb64)    (byte)limb64,(byte)(limb64>>8),\
                            (byte)(limb64>>16),(byte)(limb64>>24),\
                            (byte)(limb64>>32),(byte)(limb64>>40),\
                            (byte)(limb64>>48),(byte)(limb64>>56)
typedef byte pow256[256 / 8];

/*
 * Internal Boolean type, Boolean by value, hence safe to cast to or
 * reinterpret as 'bool'.
 */
typedef limb_t bool_t;



#ifdef __UINTPTR_TYPE__
typedef __UINTPTR_TYPE__ uptr_t;
#else
typedef const void *uptr_t;
#endif

#if !defined(restrict)
# if !defined(__STDC_VERSION__) || __STDC_VERSION__<199901
#  if defined(__GNUC__) && __GNUC__>=2
#   define restrict __restrict__
#  elif defined(_MSC_VER)
#   define restrict __restrict
#  else
#   define restrict
#  endif
# endif
#endif

#if !defined(inline) && !defined(__cplusplus)
# if !defined(__STDC_VERSION__) || __STDC_VERSION__<199901
#  if defined(__GNUC__) && __GNUC__>=2
#   define inline __inline__
#  elif defined(_MSC_VER)
#   define inline __inline
#  else
#   define inline
#  endif
# endif
#endif

#if defined(__GNUC__) || defined(__clang__)
# define BLST_launder(var) __asm__ __volatile__("" : "+r"(var))
#else
# define BLST_launder(var)
#endif

#endif

