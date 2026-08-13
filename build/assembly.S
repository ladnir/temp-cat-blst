/*
 * Copyright Supranational LLC
 * Copyright 2026 Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * Aggregate selector for the checked-in assembly providers. Platform and
 * capability macros are public build inputs; they do not depend on secrets.
 */

#include "config.h"

#ifdef BLST_PORTABLE
#define __BLST_PORTABLE__
#endif

#if defined(__x86_64) || defined(__x86_64__)
# if defined(__ELF__) 
#  include "elf/sha256-x86_64.s"
#
#  if defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "elf/ctx_inverse_mod_384-x86_64.s"
#  endif
#  if !defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "elf/ctq_inverse_mod_384-x86_64.s"
#  endif
#
#  include "elf/add_mod_384-x86_64.s"
#  include "elf/add_mod_384x384-x86_64.s"
#
#  if defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "elf/mulx_mont_384-x86_64.s"
#   include "elf/mulx_mont_256-x86_64.s"
#  endif
#  if !defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "elf/mulq_mont_384-x86_64.s"
#   include "elf/mulq_mont_256-x86_64.s"
#  endif
#
#  include "elf/add_mod_256-x86_64.s"
#  include "elf/ct_inverse_mod_256-x86_64.s"
#  include "elf/div3w-x86_64.s"
#  include "elf/ct_is_square_mod_384-x86_64.s"
#
#
# elif defined(_WIN64) || defined(__CYGWIN__)
#  include "coff/sha256-x86_64.s"
#  if defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "coff/ctx_inverse_mod_384-x86_64.s"
#  endif
#  if !defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "coff/ctq_inverse_mod_384-x86_64.s"
#  endif
#  include "coff/add_mod_384-x86_64.s"
#  include "coff/add_mod_384x384-x86_64.s"

#  if defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "coff/mulx_mont_384-x86_64.s"
#   include "coff/mulx_mont_256-x86_64.s"
#  endif
#  if !defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "coff/mulq_mont_384-x86_64.s"
#   include "coff/mulq_mont_256-x86_64.s"
#  endif
#  include "coff/add_mod_256-x86_64.s"
#  include "coff/ct_inverse_mod_256-x86_64.s"
#  include "coff/div3w-x86_64.s"
#  include "coff/ct_is_square_mod_384-x86_64.s"
# elif defined(__APPLE__)
#  include "mach-o/sha256-x86_64.s"
#  if defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "mach-o/ctx_inverse_mod_384-x86_64.s"
#  endif
#  if !defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "mach-o/ctq_inverse_mod_384-x86_64.s"
#  endif
#  include "mach-o/add_mod_384-x86_64.s"
#  include "mach-o/add_mod_384x384-x86_64.s"
#  if defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "mach-o/mulx_mont_384-x86_64.s"
#   include "mach-o/mulx_mont_256-x86_64.s"
#  endif
#  if !defined(BLST_ADX) || defined(__BLST_PORTABLE__)
#   include "mach-o/mulq_mont_384-x86_64.s"
#   include "mach-o/mulq_mont_256-x86_64.s"
#  endif
#  include "mach-o/add_mod_256-x86_64.s"
#  include "mach-o/ct_inverse_mod_256-x86_64.s"
#  include "mach-o/div3w-x86_64.s"
#  include "mach-o/ct_is_square_mod_384-x86_64.s"
# endif
#elif defined(__aarch64__)
# if defined(__CHERI_PURE_CAPABILITY__)
#  include "cheri/sha256-armv8.S"
#  include "cheri/ct_inverse_mod_384-armv8.S"
#  include "cheri/add_mod_384-armv8.S"
#  define __add_mod_384     __add_mont_384
#  define __sub_mod_384     __sub_mont_384
#  include "cheri/mul_mont_384-armv8.S"
#  include "cheri/mul_mont_256-armv8.S"
#  include "cheri/add_mod_256-armv8.S"
#  include "cheri/ct_inverse_mod_256-armv8.S"
#  include "cheri/div3w-armv8.S"
#  include "cheri/ct_is_square_mod_384-armv8.S"
# elif defined(__ELF__)
#  include "elf/sha256-armv8.S"
#  include "elf/ct_inverse_mod_384-armv8.S"
#  include "elf/add_mod_384-armv8.S"
#  define __add_mod_384     __add_mont_384
#  define __sub_mod_384     __sub_mont_384
#  include "elf/mul_mont_384-armv8.S"
#  include "elf/mul_mont_256-armv8.S"
#  include "elf/add_mod_256-armv8.S"
#  include "elf/ct_inverse_mod_256-armv8.S"
#  include "elf/div3w-armv8.S"
#  include "elf/ct_is_square_mod_384-armv8.S"
# elif defined(_WIN64)
#  include "coff/sha256-armv8.S"
#  include "coff/ct_inverse_mod_384-armv8.S"
#  include "coff/add_mod_384-armv8.S"
#  define __add_mod_384     __add_mont_384
#  define __sub_mod_384     __sub_mont_384
#  include "coff/mul_mont_384-armv8.S"
#  include "coff/mul_mont_256-armv8.S"
#  include "coff/add_mod_256-armv8.S"
#  include "coff/ct_inverse_mod_256-armv8.S"
#  include "coff/div3w-armv8.S"
#  include "coff/ct_is_square_mod_384-armv8.S"
# elif defined(__APPLE__)
#  include "mach-o/sha256-armv8.S"
#  include "mach-o/ct_inverse_mod_384-armv8.S"
#  include "mach-o/add_mod_384-armv8.S"
#  define __add_mod_384     __add_mont_384
#  define __sub_mod_384     __sub_mont_384
#  include "mach-o/mul_mont_384-armv8.S"
#  include "mach-o/mul_mont_256-armv8.S"
#  include "mach-o/add_mod_256-armv8.S"
#  include "mach-o/ct_inverse_mod_256-armv8.S"
#  include "mach-o/div3w-armv8.S"
#  include "mach-o/ct_is_square_mod_384-armv8.S"
# endif
#elif !defined(BLST_ASM)
/* The build disabled assembly; this translation unit intentionally is empty. */
#else
# error "unsupported assembly platform: Category-blst requires a supported 64-bit target"
#endif
