#ifndef __BLS12_381_DEFINES_H__
#define __BLS12_381_DEFINES_H__

#include "config.h"

// uint8_t, etc
#include <stdint.h>

#if !defined(BLST_restrict)
# if !defined(__STDC_VERSION__) || __STDC_VERSION__<199901
#  if defined(__GNUC__) && __GNUC__>=2
#   define BLST_restrict __restrict__
#  elif defined(_MSC_VER)
#   define BLST_restrict __restrict
#  else
#   define BLST_restrict
#  endif
# endif
#endif


// OC_FORCEINLINE ---------------------------------------------//
// Macro to use in place of 'inline' to force a function to be inline
#if !defined(BLST_FORCEINLINE)
#  if defined(_MSC_VER)
#    define BLST_FORCEINLINE __forceinline
#  elif defined(__GNUC__) && __GNUC__ > 3
#    define BLST_FORCEINLINE inline __attribute__ ((__always_inline__))
#  else
#    define BLST_FORCEINLINE static inline
#  endif
#endif

#  if defined(__CUDA_ARCH__) && defined(BLST_FAST_BUILD)
// cuda inlines aggretively and its results in slow build. 
#define BLST_INLINE  inline __noinline__
#else
#define BLST_INLINE  inline
#endif

#if defined(__GNUC__) || defined(__clang__)
# define BLST_launder(var) __asm__ __volatile__("" : "+r"(var))
#elif defined(_MSC_VER) && !defined(__CUDA_ARCH__)
# include <intrin.h>
# define BLST_launder(var) do { _ReadWriteBarrier(); (void)(var); } while (0)
#else
# define BLST_launder(var)
#endif

#if !defined(BLST_ASSERT)
#  include <assert.h>
#  define BLST_ASSERT(condition) assert(condition)
# endif



# if defined(__GNUC__) && __GNUC__>=5
#  define BLST_ASSUME(condition) do{BLST_ASSERT(condition); if (!(condition)) __builtin_unreachable();} while(0)
# elif defined(_MSC_VER)
#  define BLST_ASSUME(condition) do { BLST_ASSERT(condition); __assume(condition) ;} while(0)
# else
#  define BLST_ASSUME(condition) BLST_ASSERT(condition)
# endif


/*
 * Some compilers get arguably overzealous(*) when passing pointer to
 * multi-dimensional array [such as Fp2] as 'const' argument.
 * General direction seems to be to legitimize such constification,
 * so it's argued that suppressing the warning is appropriate.
 *
 * (*)  http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1923.htm
 */
#if defined(__INTEL_COMPILER)
# pragma warning(disable:167)
# pragma warning(disable:556)
#elif defined(__GNUC__) && !defined(__clang__) && (__STDC_VERSION__-0) < 202311
# pragma GCC diagnostic ignored "-Wpedantic"
#elif defined(_MSC_VER)
# pragma warning(disable: 4127 4189)
#endif

#if !defined(__wasm__) && __STDC_HOSTED__-0 != 0
# include <stdlib.h>
#endif

#if defined(__GNUC__)
# ifndef alloca
#  define alloca(s) __builtin_alloca(s)
# endif
#elif defined(__sun)
# include <alloca.h>
#elif defined(_WIN32)
# include <malloc.h>
# ifndef alloca
#  define alloca(s) _alloca(s)
# endif
#endif

#if BLST_LIMB_SIZE == 64
# define TO_LIMB_T(limb64)     limb64
#else
# define TO_LIMB_T(limb64)     (limb_t)limb64,(limb_t)(limb64>>32)
#endif




namespace blst
{

#if BLST_LIMB_SIZE == 64
    typedef unsigned long long limb_t;
#elif BLST_LIMB_SIZE == 32
    typedef unsigned int limb_t;
#else
# error "Unsupported BLST_LIMB_SIZE"
#endif
    typedef unsigned char byte;


	using u8 = uint8_t;
	using u32 = uint32_t;
	using u64 = uint64_t;
	
	using i8 = int8_t;
	using i32 = int32_t;
	using i64 = int64_t;
    typedef limb_t bool_t;

    constexpr int NLIMBS(int bits) { return bits / BLST_LIMB_SIZE; }

    typedef limb_t vec_256[NLIMBS(256)];
    typedef limb_t vec_384[NLIMBS(384)];
    typedef limb_t vec_512[NLIMBS(512)];
    typedef limb_t vec_768[NLIMBS(768)];
    typedef vec_384 vec_384x[2];      /* 0 is "real" part, 1 is "imaginary" */
    typedef vec_768 vec_768x[2];      /* 0 is "real" part, 1 is "imaginary" */

    BLST_HD void vec_zero(void* ret, size_t num);



    BLST_HD bool_t vec_is_equal(const void* a, const void* b, size_t num);
	// 512 bits, unreduced "Fr element"
    struct UFr
    {
        vec_512 data;

        BLST_HD operator limb_t* () { return data; }
        BLST_HD operator const limb_t* () const { return data; }
    };

	// 768 bits, unreduced "Fp element"
    struct UFp
    {
        vec_768 data;

		BLST_HD operator limb_t* () { return data; }
        BLST_HD operator const limb_t* () const { return data; }
    };

	// 2 * 768 bits, unreduced "Fp2 element"
    struct UFp2
    {
        UFp data[2];

        BLST_HD UFp& operator[](size_t idx) { return data[idx]; }
        BLST_HD const UFp& operator[](size_t idx) const { return data[idx]; }

    };

    //using vec_256 = Fr&;
    //using vec_384 = Fp&;
    //using vec_384x = Fp2&;
    //using vec_512 = UFr&;
    //using vec_768 = UFp&;
    //using vec_768x = UFpx2&;



#define TO_BYTES(limb64)    (byte)limb64,(byte)(limb64>>8),\
                            (byte)(limb64>>16),(byte)(limb64>>24),\
                            (byte)(limb64>>32),(byte)(limb64>>40),\
                            (byte)(limb64>>48),(byte)(limb64>>56)
    typedef byte pow256[256 / 8];


#ifdef __UINTPTR_TYPE__
    typedef __UINTPTR_TYPE__ uptr_t;
#else
    typedef const void* uptr_t;
#endif




#if defined(__INTELLISENSE__) && defined(BLST_CUDA)
 // IntelliSense: pretend we're compiling for CUDA
    struct dim3 {
        unsigned int x, y, z;
        dim3() {}
        dim3(auto, auto, auto) {}
    } ;
    static const dim3 blockIdx, threadIdx, blockDim;
#endif


}

#endif

