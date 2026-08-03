/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * This CatBLST header is locally composed. Its representation and compiler
 * support lineage includes the upstream BLST src/vect.h implementation.
 */

#ifndef __BLS12_381_DEFINES_H__
#define __BLS12_381_DEFINES_H__

#include "config.h"

#if !defined(BLST_DEFAULT_AVX)
# if defined(BLST_AVX512)
#  define BLST_DEFAULT_AVX 1
# else
#  define BLST_DEFAULT_AVX 0
# endif
#endif

// uint8_t, etc
#include <stdint.h>
#include <bit>

// For custom assert implementation
#ifndef NDEBUG
#  include <cstdio>  // for printf
#  include <cstdlib> // for abort
#  if defined(_MSC_VER) && !defined(__CUDA_ARCH__)
extern "C" __declspec(dllimport) int __stdcall IsDebuggerPresent(void);
#  endif
#endif

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
#    define BLST_FORCEINLINE inline
#  endif
#endif

#  if defined(__CUDA_ARCH__) && defined(BLST_FAST_BUILD)
// cuda inlines aggretively and its results in slow build. 
#define BLST_INLINE  inline __noinline__
#else
#define BLST_INLINE  inline
#endif

// used for unity to add the inline 
#ifndef BLST_CPP_INLINE
#define BLST_CPP_INLINE
#endif

#if defined(__GNUC__) || defined(__clang__)
# if defined(__NVCC__) || defined(__CUDACC__)
#  define BLST_launder(var) do { (void)(var); } while(0)
# else
#  define BLST_launder(var) __asm__ __volatile__("" : "+r"(var))
# endif
#elif defined(_MSC_VER) && !defined(__CUDA_ARCH__)
# include <intrin.h>
# define BLST_launder(var) do { _ReadWriteBarrier(); (void)(var); } while (0)
#else
# define BLST_launder(var)
#endif

#if !defined(BLST_ASSERT)
#  ifdef NDEBUG
     // In release builds, disable assertions
#    define BLST_ASSERT(condition) ((void)0)
#  else
     // Custom assert implementation for debug builds
namespace blst {
    namespace detail {
#if defined(__CUDA_ARCH__)
        // Device code version - simpler implementation
        __device__ BLST_FORCEINLINE void device_assert_handler(const char* condition_str, const char* file, int line) {
            printf("BLST_ASSERT failed on GPU: %s at %s:%d\n", condition_str, file, line);
            __trap(); // CUDA device trap - will break into debugger if attached
        }
#else
        // Host code version with full debugging support
        BLST_FORCEINLINE void assert_handler(const char* condition_str, const char* file, int line, const char* function) {
            // Print detailed assertion information
            fprintf(stderr, "BLST_ASSERT failed: %s\n", condition_str);
            fprintf(stderr, "  File: %s\n", file);
            fprintf(stderr, "  Line: %d\n", line);
            fprintf(stderr, "  Function: %s\n", function);
            fflush(stderr);

#if defined(_MSC_VER)
            // Visual Studio specific - trigger debugger breakpoint
            if (IsDebuggerPresent()) {
                __debugbreak();
            }
            else {
                abort();
            }
#elif defined(__GNUC__) || defined(__clang__)
            // GCC/Clang - use builtin trap for debugger
            __builtin_trap();
#else
            // Fallback - standard abort
            abort();
#endif
        }
#endif
    }
}

#if defined(__CUDA_ARCH__)
// Device code macro
#    define BLST_ASSERT(condition) \
       do { \
           if (!(condition)) { \
               ::blst::detail::device_assert_handler(#condition, __FILE__, __LINE__); \
           } \
       } while(0)
#else
// Host code macro
#    if defined(__GNUC__) || defined(__clang__)
#      define BLST_ASSERT(condition) \
         do { \
             if (!(condition)) { \
                 ::blst::detail::assert_handler(#condition, __FILE__, __LINE__, __PRETTY_FUNCTION__); \
             } \
         } while(0)
#    elif defined(_MSC_VER)
#      define BLST_ASSERT(condition) \
         do { \
             if (!(condition)) { \
                 ::blst::detail::assert_handler(#condition, __FILE__, __LINE__, __FUNCTION__); \
             } \
         } while(0)
#    else
#      define BLST_ASSERT(condition) \
         do { \
             if (!(condition)) { \
                 ::blst::detail::assert_handler(#condition, __FILE__, __LINE__, "unknown"); \
             } \
         } while(0)
#    endif
#endif
#  endif
#endif


# if defined(__GNUC__) && __GNUC__>=5
#  define BLST_ASSUME(condition) do{BLST_ASSERT(condition); if (!(condition)) __builtin_unreachable();} while(0)
# elif defined(_MSC_VER)
#  define BLST_ASSUME(condition) do { BLST_ASSERT(condition); __assume(condition) ;} while(0)
# else
#  define BLST_ASSUME(condition) BLST_ASSERT(condition)
# endif


#if !defined(__wasm__) && __STDC_HOSTED__-0 != 0
# include <stdlib.h>
#endif

#if BLST_LIMB_SIZE != 64
# error "CatBLST supports only 64-bit limbs"
#endif
#define TO_LIMB_T(limb64)     (limb_t)limb64


#if defined(__CUDACC__)
 /* CUDA build: make these no-ops to avoid pragma issues in shared headers. */
#define BLST_UNROLL        /* no-op */
#elif defined(_MSC_VER)

 /* MSVC: use __pragma */
#define BLST_UNROLL        __pragma(loop(ivdep))
#else

 /* GCC/Clang: use _Pragma (and keep it generic). */
#define BLST_UNROLL        _Pragma(unroll)
#endif


namespace blst
{

    using limb_t = uint64_t;
    typedef unsigned char byte;


	using u8 = uint8_t;
	using u16 = uint16_t;
	using u32 = uint32_t;
	using u64 = uint64_t;
	
	using i8 = int8_t;
	using i16 = int16_t;
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

    inline constexpr u64 log2ceil(u64 x)
    {
        BLST_ASSERT(x != 0);
        if (x == 0) {
            return 0;
        }
        return std::bit_width(x - 1);
    }


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

#define TO_BYTES(limb64)    (byte)limb64,(byte)(limb64>>8),\
                            (byte)(limb64>>16),(byte)(limb64>>24),\
                            (byte)(limb64>>32),(byte)(limb64>>40),\
                            (byte)(limb64>>48),(byte)(limb64>>56)
    typedef byte pow256[256 / 8];


    typedef uintptr_t uptr_t;


    // Zero is the fail-closed Release result for an invalid zero divisor.
    inline constexpr size_t ceil_div(size_t x, size_t y)
    {
        BLST_ASSERT(y != 0);
        return y == 0 ? 0 : x / y + (x % y != 0);
    }


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

