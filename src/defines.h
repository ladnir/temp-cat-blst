#ifndef __BLS12_381_DEFINES_H__
#define __BLS12_381_DEFINES_H__

#include "config.h"

#if defined(__x86_64__) || defined(__aarch64__)
/* These are available even in ILP32 flavours, but even then they are
 * capable of performing 64-bit operations as efficiently as in *P64. */
typedef unsigned long long limb_t;
# define LIMB_T_BITS    64

#elif defined(_WIN64)   /* Win64 is P64 */

#ifdef __BLST_NO_ASM__

typedef unsigned __int32 limb_t;
# define LIMB_T_BITS    32
#else
typedef unsigned __int64 limb_t;
# define LIMB_T_BITS    64
#endif // __BLST_NO_ASM__

#elif defined(__BLST_NO_ASM__) || defined(__wasm64__)
typedef unsigned int limb_t;
# define LIMB_T_BITS    32
# ifndef __BLST_NO_ASM__
#  define __BLST_NO_ASM__
# endif

#else                   /* 32 bits on 32-bit platforms, 64 - on 64-bit */
typedef unsigned long limb_t;
#  ifdef _LP64
#   define LIMB_T_BITS   64
#  else
#   define LIMB_T_BITS   32
#   define __BLST_NO_ASM__
#  endif
#endif
typedef unsigned char byte;
/*
 * Why isn't LIMB_T_BITS defined as 8*sizeof(limb_t)? Because pre-processor
 * knows nothing about sizeof(anything)...
 */
#if LIMB_T_BITS == 64
# define TO_LIMB_T(limb64)     limb64
#else
# define TO_LIMB_T(limb64)     (limb_t)limb64,(limb_t)(limb64>>32)
#endif

#define NLIMBS(bits)   (bits/LIMB_T_BITS)

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

