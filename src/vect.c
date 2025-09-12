/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "vect.h"
 //
 //#ifndef BLST_ASM
 //# include "no_asm.h"
 //#endif

 /*
  * Following are some reference C implementations to assist new
  * assembly modules development, as starting-point stand-ins and for
  * cross-checking. In order to "polyfil" specific subroutine redefine
  * it on compiler command line, e.g. -Dmul_mont_384x=_mul_mont_384x.
  */



  //#ifndef BLST_ASM
  //#endif

  //#ifndef BLST_ASM
  //inline void mul_by_3_mod_384(vec384 ret, const vec384 a, const vec384 mod)
  //{
  //    vec384 t;
  //
  //    add_mod_384(t, a, a, mod);
  //    add_mod_384(ret, t, a, mod);
  //}
  //#endif

  //
  //#if defined(__BLST_NO_ASM__) && !(defined(__ADX__) && !defined(__BLST_PORTABLE__))
  //void mul_mont_384x(vec384x ret, const vec384x a, const vec384x b,
  //                   const vec384 mod, limb_t n0)
  //{
  //    vec768 t0, t1, t2;
  //    vec384 aa, bb;
  //
  //    mul_384(t0, a[0], b[0]);
  //    mul_384(t1, a[1], b[1]);
  //
  //    add_mod_384(aa, a[0], a[1], mod);
  //    add_mod_384(bb, b[0], b[1], mod);
  //    mul_384(t2, aa, bb);
  //    sub_mod_384x384(t2, t2, t0, mod);
  //    sub_mod_384x384(t2, t2, t1, mod);
  //
  //    sub_mod_384x384(t0, t0, t1, mod);
  //
  //    redc_mont_384(ret[0], t0, mod, n0);
  //    redc_mont_384(ret[1], t2, mod, n0);
  //}
  //#endif
  //

  //void mul_mont_384x(vec384x ret, const vec384x a, const vec384x b,
  //    const vec384 p, limb_t n0)
  //{
  //    vec384 aa, bb, cc;
  //
  //    add_mod_n(aa, a[0], a[1], p, NLIMBS(384));
  //    add_mod_n(bb, b[0], b[1], p, NLIMBS(384));
  //    mul_mont_n(bb, bb, aa, p, n0, NLIMBS(384));
  //    mul_mont_n(aa, a[0], b[0], p, n0, NLIMBS(384));
  //    mul_mont_n(cc, a[1], b[1], p, n0, NLIMBS(384));
  //    sub_mod_n(ret[0], aa, cc, p, NLIMBS(384));
  //    sub_mod_n(ret[1], bb, aa, p, NLIMBS(384));
  //    sub_mod_n(ret[1], ret[1], cc, p, NLIMBS(384));
  //}

