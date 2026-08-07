#pragma once

#if defined(_MSC_VER) && !defined(__clang__)
#define INLINE __forceinline
#define NOINLINE __declspec(noinline)
#else
#define INLINE inline __attribute__((always_inline))
#define NOINLINE __attribute__((noinline))
#endif

// Flags controlling explicit options
// Set to INLINE or NOINLINE for explicit control
#define EXPAND_IS_INLINE INLINE
#define REDUCE_IS_INLINE INLINE
// Used to see code in compiler explorer
// #define TEST_INLINE inline
#define TEST_INLINE INLINE
