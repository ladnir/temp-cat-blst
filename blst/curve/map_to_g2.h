/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_MAP_TO_G2_H
#define BLST_MAP_TO_G2_H

#include "blst/curve/point.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/hash_to_field.h"

namespace blst
{

    //BLST_EXTERN_C_BEGIN
    // /*
    //  * y^2 = x^3 + A'*x + B', isogenous one
    //  */
    //extern const Fp2 Aprime_E2;
    //extern const Fp2 Bprime_E2;

    //BLST_HD void map_fp2_times_Zz(Fp2 map[], const Fp2 isogeny_map[],
    //    const Fp2 Zz_powers[], size_t n);
    //BLST_HD void map_fp2(Fp2& acc, const Fp2& x, const Fp2& map, size_t n);
    //BLST_HD void isogeny_map_to_E2(POINTonE2& out, const POINTonE2& p);
    //BLST_HD void map_to_isogenous_E2(POINTonE2& p, const Fp2& u);

    ///*
    // * As per suggestions in "7. Clearing the cofactor" at
    // * https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06
    // */
    //BLST_HD void POINTonE2_add_n_dbl(POINTonE2& out, const POINTonE2& p, size_t n);

    //BLST_HD void POINTonE2_times_minus_z(POINTonE2& out, const POINTonE2& in);

    ////void psi(POINTonE2& out, const POINTonE2& in);

    //BLST_HD void clear_cofactor(POINTonE2& out, const POINTonE2& p);
    ///*
    // * |u|, |v| are expected to be in Montgomery representation
    // */
    //BLST_HD void map_to_g2(POINTonE2& out, const Fp2& u, const Fp2& v);

    //BLST_HD void blst_map_to_g2(POINTonE2& out, const Fp2& u, const Fp2& v);

    //BLST_HD void Encode_to_G2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
    //    const unsigned char* DST, size_t DST_len,
    //    const unsigned char* aug, size_t aug_len);

    //BLST_HD void blst_encode_to_g2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
    //    const unsigned char* DST, size_t DST_len,
    //    const unsigned char* aug, size_t aug_len);

    //BLST_HD void Hash_to_G2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
    //    const unsigned char* DST, size_t DST_len,
    //    const unsigned char* aug, size_t aug_len);

    //BLST_HD void blst_hash_to_g2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
    //    const unsigned char* DST, size_t DST_len,
    //    const unsigned char* aug, size_t aug_len);

    //BLST_HD bool_t POINTonE2_in_G2(const POINTonE2& P);

    //BLST_HD int blst_p2_in_g2(const POINTonE2& p);

    //BLST_HD int blst_p2_affine_in_g2(const POINTonE2_affine* p);
















    /*
     * y^2 = x^3 + A'*x + B', isogenous one
     */
#define BLST_Aprime_E2_LIMBS      /* 240*i */                       \
  Fp{ 0 },                                                            \
  Fp{ TO_LIMB_T(0xe53a000003135242), TO_LIMB_T(0x01080c0fdef80285),   \
    TO_LIMB_T(0xe7889edbe340f6bd), TO_LIMB_T(0x0b51375126310601),   \
    TO_LIMB_T(0x02d6985717c744ab), TO_LIMB_T(0x1220b4e979ea5467) }

#define BLST_Bprime_E2LIMBS      /* 1012 + 1012*i */                \
  Fp{ TO_LIMB_T(0x22ea00000cf89db2), TO_LIMB_T(0x6ec832df71380aa4),   \
    TO_LIMB_T(0x6e1b94403db5a66e), TO_LIMB_T(0x75bf3c53a79473ba),   \
    TO_LIMB_T(0x3dd3a569412c0a34), TO_LIMB_T(0x125cdb5e74dc4fd1) }, \
  Fp{ TO_LIMB_T(0x22ea00000cf89db2), TO_LIMB_T(0x6ec832df71380aa4),   \
    TO_LIMB_T(0x6e1b94403db5a66e), TO_LIMB_T(0x75bf3c53a79473ba),   \
    TO_LIMB_T(0x3dd3a569412c0a34), TO_LIMB_T(0x125cdb5e74dc4fd1) }

     /*
      * y^2 = x^3 + A'*x + B', isogenous one
      */
    static const Fp2 BLST_Aprime_E2_HOST = { BLST_Aprime_E2_LIMBS };
    static const Fp2 BLST_Bprime_E2_HOST = { BLST_Bprime_E2LIMBS };

#if defined(BLST_CUDA) && defined(__CUDACC__)

    /*
     * y^2 = x^3 + A'*x + B', isogenous one
     */
    static __constant__ Fp2 BLST_Aprime_E2_DEV = { BLST_Aprime_E2_LIMBS };
    static __constant__ Fp2 BLST_Bprime_E2_DEV = { BLST_Bprime_E2LIMBS };

#endif


    // the A' in  y^2 = x ^ 3 + A'*x + B', isogenous one
    BLST_HD BLST_FORCEINLINE const Fp2& get_BLST_Aprime_E2(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_Aprime_E2_DEV;
#else
        return BLST_Aprime_E2_HOST;
#endif
    }

    // the B' in  y^2 = x ^ 3 + A'*x + B', isogenous one
    BLST_HD BLST_FORCEINLINE const Fp2& get_BLST_Bprime_E2(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_Bprime_E2_DEV;
#else
        return BLST_Bprime_E2_HOST;
#endif
    }









    BLST_HD BLST_INLINE void map_fp2_times_Zz(Fp2 map[], const Fp2 isogeny_map[],
        const Fp2 Zz_powers[], size_t n)
    {
        while (n--)
            Fp2::mul(map[n], isogeny_map[n], Zz_powers[n]);
    }

    BLST_HD BLST_INLINE void map_fp2(Fp2& acc, const Fp2& x, const Fp2 map[], size_t n)
    {
        while (n--) {
            Fp2::mul(acc, acc, x);
            Fp2::add(acc, acc, map[n]);
        }
    }

    BLST_HD BLST_INLINE void isogeny_map_to_E2(POINTonE2& out, const POINTonE2& p)
    {
        /*
         * x = x_num / x_den, where
         * x_num = k_(1,3) * x'^3 + k_(1,2) * x'^2 + k_(1,1) * x' + k_(1,0)
         * ...
         */
        const Fp2 isogeny_map_x_num[] = {    /* (k_(1,*)<<384) % P   */
         {Fp{ TO_LIMB_T(0x47f671c71ce05e62), TO_LIMB_T(0x06dd57071206393e),
            TO_LIMB_T(0x7c80cd2af3fd71a2), TO_LIMB_T(0x048103ea9e6cd062),
            TO_LIMB_T(0xc54516acc8d037f6), TO_LIMB_T(0x13808f550920ea41) },
          Fp{ TO_LIMB_T(0x47f671c71ce05e62), TO_LIMB_T(0x06dd57071206393e),
            TO_LIMB_T(0x7c80cd2af3fd71a2), TO_LIMB_T(0x048103ea9e6cd062),
            TO_LIMB_T(0xc54516acc8d037f6), TO_LIMB_T(0x13808f550920ea41) }},
         {Fp{ 0 },
          Fp{ TO_LIMB_T(0x5fe55555554c71d0), TO_LIMB_T(0x873fffdd236aaaa3),
            TO_LIMB_T(0x6a6b4619b26ef918), TO_LIMB_T(0x21c2888408874945),
            TO_LIMB_T(0x2836cda7028cabc5), TO_LIMB_T(0x0ac73310a7fd5abd) }},
         {Fp{ TO_LIMB_T(0x0a0c5555555971c3), TO_LIMB_T(0xdb0c00101f9eaaae),
            TO_LIMB_T(0xb1fb2f941d797997), TO_LIMB_T(0xd3960742ef416e1c),
            TO_LIMB_T(0xb70040e2c20556f4), TO_LIMB_T(0x149d7861e581393b) },
          Fp{ TO_LIMB_T(0xaff2aaaaaaa638e8), TO_LIMB_T(0x439fffee91b55551),
            TO_LIMB_T(0xb535a30cd9377c8c), TO_LIMB_T(0x90e144420443a4a2),
            TO_LIMB_T(0x941b66d3814655e2), TO_LIMB_T(0x0563998853fead5e) }},
         {Fp{ TO_LIMB_T(0x40aac71c71c725ed), TO_LIMB_T(0x190955557a84e38e),
            TO_LIMB_T(0xd817050a8f41abc3), TO_LIMB_T(0xd86485d4c87f6fb1),
            TO_LIMB_T(0x696eb479f885d059), TO_LIMB_T(0x198e1a74328002d2) },
          Fp{ 0 }}
        };
        /* ...
         * x_den = x'^2 + k_(2,1) * x' + k_(2,0)
         */
        const Fp2 isogeny_map_x_den[] = {    /* (k_(2,*)<<384) % P   */
         {Fp{ 0 },
          Fp{ TO_LIMB_T(0x1f3affffff13ab97), TO_LIMB_T(0xf25bfc611da3ff3e),
            TO_LIMB_T(0xca3757cb3819b208), TO_LIMB_T(0x3e6427366f8cec18),
            TO_LIMB_T(0x03977bc86095b089), TO_LIMB_T(0x04f69db13f39a952) }},
         {Fp{ TO_LIMB_T(0x447600000027552e), TO_LIMB_T(0xdcb8009a43480020),
            TO_LIMB_T(0x6f7ee9ce4a6e8b59), TO_LIMB_T(0xb10330b7c0a95bc6),
            TO_LIMB_T(0x6140b1fcfb1e54b7), TO_LIMB_T(0x0381be097f0bb4e1) },
          Fp{ TO_LIMB_T(0x7588ffffffd8557d), TO_LIMB_T(0x41f3ff646e0bffdf),
            TO_LIMB_T(0xf7b1e8d2ac426aca), TO_LIMB_T(0xb3741acd32dbb6f8),
            TO_LIMB_T(0xe9daf5b9482d581f), TO_LIMB_T(0x167f53e0ba7431b8) }}
        };
        /*
         * y = y' * y_num / y_den, where
         * y_num = k_(3,3) * x'^3 + k_(3,2) * x'^2 + k_(3,1) * x' + k_(3,0)
         * ...
         */
        const Fp2 isogeny_map_y_num[] = {    /* (k_(3,*)<<384) % P   */
         {Fp{ TO_LIMB_T(0x96d8f684bdfc77be), TO_LIMB_T(0xb530e4f43b66d0e2),
            TO_LIMB_T(0x184a88ff379652fd), TO_LIMB_T(0x57cb23ecfae804e1),
            TO_LIMB_T(0x0fd2e39eada3eba9), TO_LIMB_T(0x08c8055e31c5d5c3) },
          Fp{ TO_LIMB_T(0x96d8f684bdfc77be), TO_LIMB_T(0xb530e4f43b66d0e2),
            TO_LIMB_T(0x184a88ff379652fd), TO_LIMB_T(0x57cb23ecfae804e1),
            TO_LIMB_T(0x0fd2e39eada3eba9), TO_LIMB_T(0x08c8055e31c5d5c3) }},
         {Fp{ 0 },
          Fp{ TO_LIMB_T(0xbf0a71c71c91b406), TO_LIMB_T(0x4d6d55d28b7638fd),
            TO_LIMB_T(0x9d82f98e5f205aee), TO_LIMB_T(0xa27aa27b1d1a18d5),
            TO_LIMB_T(0x02c3b2b2d2938e86), TO_LIMB_T(0x0c7d13420b09807f) }},
         {Fp{ TO_LIMB_T(0xd7f9555555531c74), TO_LIMB_T(0x21cffff748daaaa8),
            TO_LIMB_T(0x5a9ad1866c9bbe46), TO_LIMB_T(0x4870a2210221d251),
            TO_LIMB_T(0x4a0db369c0a32af1), TO_LIMB_T(0x02b1ccc429ff56af) },
          Fp{ TO_LIMB_T(0xe205aaaaaaac8e37), TO_LIMB_T(0xfcdc000768795556),
            TO_LIMB_T(0x0c96011a8a1537dd), TO_LIMB_T(0x1c06a963f163406e),
            TO_LIMB_T(0x010df44c82a881e6), TO_LIMB_T(0x174f45260f808feb) }},
         {Fp{ TO_LIMB_T(0xa470bda12f67f35c), TO_LIMB_T(0xc0fe38e23327b425),
            TO_LIMB_T(0xc9d3d0f2c6f0678d), TO_LIMB_T(0x1c55c9935b5a982e),
            TO_LIMB_T(0x27f6c0e2f0746764), TO_LIMB_T(0x117c5e6e28aa9054) },
          Fp{ 0 }}
        };
        /* ...
         * y_den = x'^3 + k_(4,2) * x'^2 + k_(4,1) * x' + k_(4,0)
         */
        const Fp2 isogeny_map_y_den[] = {    /* (k_(4,*)<<384) % P   */
         {Fp{ TO_LIMB_T(0x0162fffffa765adf), TO_LIMB_T(0x8f7bea480083fb75),
            TO_LIMB_T(0x561b3c2259e93611), TO_LIMB_T(0x11e19fc1a9c875d5),
            TO_LIMB_T(0xca713efc00367660), TO_LIMB_T(0x03c6a03d41da1151) },
          Fp{ TO_LIMB_T(0x0162fffffa765adf), TO_LIMB_T(0x8f7bea480083fb75),
            TO_LIMB_T(0x561b3c2259e93611), TO_LIMB_T(0x11e19fc1a9c875d5),
            TO_LIMB_T(0xca713efc00367660), TO_LIMB_T(0x03c6a03d41da1151) }},
         {Fp{ 0 },
          Fp{ TO_LIMB_T(0x5db0fffffd3b02c5), TO_LIMB_T(0xd713f52358ebfdba),
            TO_LIMB_T(0x5ea60761a84d161a), TO_LIMB_T(0xbb2c75a34ea6c44a),
            TO_LIMB_T(0x0ac6735921c1119b), TO_LIMB_T(0x0ee3d913bdacfbf6) }},
         {Fp{ TO_LIMB_T(0x66b10000003affc5), TO_LIMB_T(0xcb1400e764ec0030),
            TO_LIMB_T(0xa73e5eb56fa5d106), TO_LIMB_T(0x8984c913a0fe09a9),
            TO_LIMB_T(0x11e10afb78ad7f13), TO_LIMB_T(0x05429d0e3e918f52) },
          Fp{ TO_LIMB_T(0x534dffffffc4aae6), TO_LIMB_T(0x5397ff174c67ffcf),
            TO_LIMB_T(0xbff273eb870b251d), TO_LIMB_T(0xdaf2827152870915),
            TO_LIMB_T(0x393a9cbaca9e2dc3), TO_LIMB_T(0x14be74dbfaee5748) }}
        };
        Fp2 Zz_powers[3], map[3], xn, xd, yn, yd;

        /* lay down Z^2 powers in descending order                          */
        Fp2::sqr(Zz_powers[2], p.Z);                       /* ZZ^1          */
        Fp2::sqr(Zz_powers[1], Zz_powers[2]);               /* ZZ^2  1+1     */
        Fp2::mul(Zz_powers[0], Zz_powers[2], Zz_powers[1]); /* ZZ^3  2+1     */

        map_fp2_times_Zz(map, isogeny_map_x_num, Zz_powers, 3);
        Fp2::mul(xn, p.X, isogeny_map_x_num[3]);
        Fp2::add(xn, xn, map[2]);
        map_fp2(xn, p.X, map, 2);

        map_fp2_times_Zz(map, isogeny_map_x_den, Zz_powers + 1, 2);
        Fp2::add(xd, p.X, map[1]);
        map_fp2(xd, p.X, map, 1);
        Fp2::mul(xd, xd, Zz_powers[2]);      /* xd *= Z^2                    */

        map_fp2_times_Zz(map, isogeny_map_y_num, Zz_powers, 3);
        Fp2::mul(yn, p.X, isogeny_map_y_num[3]);
        Fp2::add(yn, yn, map[2]);
        map_fp2(yn, p.X, map, 2);
        Fp2::mul(yn, yn, p.Y);              /* yn *= Y                      */

        map_fp2_times_Zz(map, isogeny_map_y_den, Zz_powers, 3);
        Fp2::add(yd, p.X, map[2]);
        map_fp2(yd, p.X, map, 2);
        Fp2::mul(Zz_powers[2], Zz_powers[2], p.Z);
        Fp2::mul(yd, yd, Zz_powers[2]);      /* yd *= Z^3                    */

        /* convert (xn, xd, yn, yd) to Jacobian coordinates                 */
        Fp2::mul(out.Z, xd, yd);            /* Z = xd * yd                  */
        Fp2::mul(out.X, xn, yd);
        Fp2::mul(out.X, out.X, out.Z);    /* X = xn * xd * yd^2           */
        Fp2::sqr(out.Y, out.Z);
        Fp2::mul(out.Y, out.Y, xd);
        Fp2::mul(out.Y, out.Y, yn);        /* Y = yn * xd^3 * yd^2         */
    }

    BLST_HD BLST_INLINE void map_to_isogenous_E2(POINTonE2& p, const Fp2& u)
    {
        const Fp2 minus_A = {
          Fp{ 0 },
          Fp{ TO_LIMB_T(0xd4c4fffffcec5869), TO_LIMB_T(0x1da3f3eed25bfd79),
            TO_LIMB_T(0x7fa833c5136fff67), TO_LIMB_T(0x59261433cd540cbd),
            TO_LIMB_T(0x48450f5f2b84682c), TO_LIMB_T(0x07e05d00bf959233) }
        };
        const Fp2 Z = {              /* -2 - i */
          Fp{ TO_LIMB_T(0x87ebfffffff9555c), TO_LIMB_T(0x656fffe5da8ffffa),
            TO_LIMB_T(0x0fd0749345d33ad2), TO_LIMB_T(0xd951e663066576f4),
            TO_LIMB_T(0xde291a3d41e980d3), TO_LIMB_T(0x0815664c7dfe040d) },
          Fp{ TO_LIMB_T(0x43f5fffffffcaaae), TO_LIMB_T(0x32b7fff2ed47fffd),
            TO_LIMB_T(0x07e83a49a2e99d69), TO_LIMB_T(0xeca8f3318332bb7a),
            TO_LIMB_T(0xef148d1ea0f4c069), TO_LIMB_T(0x040ab3263eff0206) }
        };
        const Fp2 ZxA = {            /* 240 - 480*i */
          Fp{ TO_LIMB_T(0xe53a000003135242), TO_LIMB_T(0x01080c0fdef80285),
            TO_LIMB_T(0xe7889edbe340f6bd), TO_LIMB_T(0x0b51375126310601),
            TO_LIMB_T(0x02d6985717c744ab), TO_LIMB_T(0x1220b4e979ea5467) },
          Fp{ TO_LIMB_T(0xa989fffff9d8b0d2), TO_LIMB_T(0x3b47e7dda4b7faf3),
            TO_LIMB_T(0xff50678a26dffece), TO_LIMB_T(0xb24c28679aa8197a),
            TO_LIMB_T(0x908a1ebe5708d058), TO_LIMB_T(0x0fc0ba017f2b2466) }
        };
        Fp2 uu, tv2, tv4, x2n, gx1, gxd, y2;
#if 0
        Fp2 xn, x1n, xd, y, y1, Zuu;
#else
# define xn     p.X
# define y      p.Y
# define xd     p.Z
# define x1n    xn
# define y1     y
# define Zuu    x2n
#endif
#define sgn0_fp2(a) (Fp2::sgn0_pty_mont(a) & 1)
        bool_t e1, e2;

        /*
         * as per map_to_curve() from poc/sswu_opt.sage at
         * https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve
         * with 9mod16 twists...
         */
         /* x numerator variants                                             */
        Fp2::sqr(uu, u);                           /* uu = u^2                     */
        Fp2::mul(Zuu, Z, uu);                      /* Zuu = Z * uu                 */
        Fp2::sqr(tv2, Zuu);                        /* tv2 = Zuu^2                  */
        Fp2::add(tv2, tv2, Zuu);                   /* tv2 = tv2 + Zuu              */
        Fp2::add(x1n, tv2, *(Fp2*)get_BLST_Rx().p2);     /* x1n = tv2 + 1                */
        Fp2::mul(x1n, x1n, get_BLST_Bprime_E2()); /* x1n = x1n * B                */
        Fp2::mul(x2n, Zuu, x1n);                   /* x2n = Zuu * x1n              */

        /* x denumenator                                                    */
        Fp2::mul(xd, minus_A, tv2);          /* xd = -A * tv2                */
        e1 = vec_is_zero(&(xd), sizeof(xd));   /* e1 = xd == 0                 */
        vec_select(&(xd), &ZxA, &(xd), sizeof(xd), e1);    /*              # If xd == 0, set xd = Z*A */

        /* y numerators variants                                            */
        Fp2::sqr(tv2, xd);                   /* tv2 = xd^2                   */
        Fp2::mul(gxd, xd, tv2);              /* gxd = xd^3                   */
        Fp2::mul(tv2, get_BLST_Aprime_E2(), tv2);       /* tv2 = A * tv2                */
        Fp2::sqr(gx1, x1n);                  /* gx1 = x1n^2                  */
        Fp2::add(gx1, gx1, tv2);             /* gx1 = gx1 + tv2      # x1n^2 + A*xd^2 */
        Fp2::mul(gx1, gx1, x1n);             /* gx1 = gx1 * x1n      # x1n^3 + A*x1n*xd^2 */
        Fp2::mul(tv2, get_BLST_Bprime_E2(), gxd);       /* tv2 = B * gxd                */
        Fp2::add(gx1, gx1, tv2);             /* gx1 = gx1 + tv2      # x1^3 + A*x1*xd^2 + B*xd^3 */
        Fp2::sqr(tv4, gxd);                  /* tv4 = gxd^2                  */
        Fp2::mul(tv2, gx1, gxd);             /* tv2 = gx1 * gxd              */
        Fp2::mul(tv4, tv4, tv2);             /* tv4 = tv4 * tv2      # gx1*gxd^3 */
        e2 = Fp2::recip_sqrt(y1, tv4);        /* y1 = tv4^c1          # (gx1*gxd^3)^((p^2-9)/16) */
        Fp2::mul(y1, y1, tv2);               /* y1 = y1 * tv2        # gx1*gxd*y1 */
        Fp2::mul(y2, y1, uu);                /* y2 = y1 * uu                 */
        Fp2::mul(y2, y2, u);                 /* y2 = y2 * u                  */

        /* choose numerators                                                */
        vec_select(&(xn), &(x1n), &x2n, sizeof(xn), e2);   /* xn = e2 ? x1n : x2n  */
        vec_select(&y, &y1, &y2, sizeof(y), e2);       /* y  = e2 ? y1 : y2    */

        e1 = sgn0_fp2(u);
        e2 = sgn0_fp2(y);
        Fp2::cneg(y, y, e1 ^ e2);              /* fix sign of y                */
        /* return (xn, xd, y, 1)        */

    /* convert (xn, xd, y, 1) to Jacobian projective coordinates        */
        Fp2::mul(p.X, xn, xd);              /* X = xn * xd                  */
        Fp2::mul(p.Y, y, gxd);              /* Y = y * xd^3                 */
#ifndef xd
        vec_copy(p.Z, xd, sizeof(xd));     /* Z = xd                       */
#else
# undef xn
# undef y
# undef xd
# undef x1n
# undef y1
# undef Zuu
# undef tv4
#endif
#undef sgn0_fp2
    }

#if 0

    BLST_HD BLST_INLINE void clear_cofactor(POINTonE2& out, const POINTonE2& p)
    {
        static const byte h_eff[] = {
            TO_BYTES(0xe8020005aaa95551), TO_BYTES(0x59894c0adebbf6b4),
            TO_BYTES(0xe954cbc06689f6a3), TO_BYTES(0x2ec0ec69d7477c1a),
            TO_BYTES(0x6d82bf015d1212b0), TO_BYTES(0x329c2f178731db95),
            TO_BYTES(0x9986ff031508ffe1), TO_BYTES(0x88e2a8e9145ad768),
            TO_BYTES(0x584c6a0ea91b3528), TO_BYTES(0x0bc69f08f2ee75b3)
        };

        POINTonE2_mult_w5(out, p, h_eff, 636);
    }
#else
    /*
     * As per suggestions in "7. Clearing the cofactor" at
     * https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06
     */
    BLST_HD BLST_INLINE void POINTonE2_add_n_dbl(POINTonE2& out, const POINTonE2& p, size_t n)
    {
        POINTonE2::dadd(out, out, p, NULL);
        while (n--)
            POINTonE2::doubling(out, out);
    }

    BLST_HD BLST_INLINE void POINTonE2_times_minus_z(POINTonE2& out, const POINTonE2& in)
    {
        POINTonE2::doubling(out, in);          /*      1: 0x2                  */
        POINTonE2_add_n_dbl(out, in, 2);    /*   2..4: 0x3..0xc             */
        POINTonE2_add_n_dbl(out, in, 3);    /*   5..8: 0xd..0x68            */
        POINTonE2_add_n_dbl(out, in, 9);    /*  9..18: 0x69..0xd200         */
        POINTonE2_add_n_dbl(out, in, 32);   /* 19..51: ..0xd20100000000     */
        POINTonE2_add_n_dbl(out, in, 16);   /* 52..68: ..0xd201000000010000 */
    }

    BLST_HD BLST_INLINE void psi(POINTonE2& out, const POINTonE2& in);

    BLST_HD BLST_INLINE void clear_cofactor(POINTonE2& out, const POINTonE2& p)
    {
        POINTonE2 t0, t1;

        /* A.Budroni, F.Pintore, "Efficient hash maps to G2 on BLS curves"  */
        POINTonE2::doubling(out, p);           /* out = 2P                     */
        psi(out, out);                      /* out = Ψ(2P)                  */
        psi(out, out);                      /* out = Ψ²(2P)                 */

        vec_copy(&t0, &p, sizeof(t0));
        POINTonE2::cneg(t0, 1);             /* t0 = -P                      */
        psi(t1, t0);                      /* t1 = -Ψ(P)                   */
        POINTonE2::dadd(out, out, t0, NULL);/* out = Ψ²(2P) - P             */
        POINTonE2::dadd(out, out, t1, NULL);/* out = Ψ²(2P) - P - Ψ(P)      */

        POINTonE2_times_minus_z(t0, p);    /* t0 = [-z]P                   */
        POINTonE2::dadd(t0, t0, p, NULL);  /* t0 = [-z + 1]P               */
        POINTonE2::dadd(t0, t0, t1, NULL);/* t0 = [-z + 1]P - Ψ(P)        */
        POINTonE2_times_minus_z(t1, t0);  /* t1 = [z² - z]P + [z]Ψ(P)     */
        POINTonE2::dadd(out, out, t1, NULL);/* out = [z² - z - 1]P          */
        /*     + [z - 1]Ψ(P)            */
        /*     + Ψ²(2P)                 */
    }
#endif

    /*
     * |u|, |v| are expected to be in Montgomery representation
     */
    BLST_HD BLST_INLINE void map_to_g2(POINTonE2& out, const Fp2& u, const Fp2* v)
    {
        POINTonE2 p;

        map_to_isogenous_E2(p, u);

        if (v != NULL) {
            map_to_isogenous_E2(out, *v);    /* borrow |out|                 */
            POINTonE2::dadd(p, p, out, &get_BLST_Aprime_E2());
        }

        isogeny_map_to_E2(p, p);          /* sprinkle isogenous powder    */
        clear_cofactor(out, p);
    }

    BLST_HD BLST_INLINE void blst_map_to_g2(POINTonE2& out, const Fp2& u, const Fp2* v)
    {
        map_to_g2(out, u, v);
    }

    BLST_HD BLST_INLINE void Encode_to_G2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len,
        const unsigned char* aug, size_t aug_len)
    {
        Fp2 u;

        hash_to_field(u.data, 2, aug, aug_len, msg, msg_len, DST, DST_len);
        map_to_g2(p, u, NULL);
    }

    BLST_HD BLST_INLINE void blst_encode_to_g2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len,
        const unsigned char* aug, size_t aug_len)
    {
        Encode_to_G2(p, msg, msg_len, DST, DST_len, aug, aug_len);
    }

    BLST_HD BLST_INLINE void Hash_to_G2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len,
        const unsigned char* aug, size_t aug_len)
    {
        Fp2 u[2];

        hash_to_field(u->data, 4, aug, aug_len, msg, msg_len, DST, DST_len);
        map_to_g2(p, u[0], &u[1]);
    }

    BLST_HD BLST_INLINE void blst_hash_to_g2(POINTonE2& p, const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len,
        const unsigned char* aug, size_t aug_len)
    {
        Hash_to_G2(p, msg, msg_len, DST, DST_len, aug, aug_len);
    }

    BLST_HD BLST_INLINE bool_t POINTonE2_in_G2(const POINTonE2& P)
    {
#if 0
        POINTonE2 t0, t1, t2;

        /* Bowe, S., "Faster subgroup checks for BLS12-381"                 */
        psi(&t0, P);                        /* Ψ(P)                         */
        psi(&t0, &t0);                      /* Ψ²(P)                        */
        psi(&t1, &t0);                      /* Ψ³(P)                        */

        POINTonE2_times_minus_z(&t2, &t1);
        POINTonE2::dadd(&t0, &t0, &t2, NULL);
        POINTonE2::cneg(&t0, 1);
        POINTonE2::dadd(&t0, &t0, P, NULL);  /* [z]Ψ³(P) - Ψ²(P) + P         */

        return vec_is_zero(t0.Z, sizeof(t0.Z));
#else
        POINTonE2 t0, t1;

        /* Scott, M., https://eprint.iacr.org/2021/1130 */
        psi(t0, P);                            /* Ψ(P) */

        POINTonE2_times_minus_z(t1, P);
        POINTonE2::cneg(t1, 1);                 /* [z]P */

        return POINTonE2::is_equal(t0, t1);
#endif
    }

    BLST_HD BLST_INLINE int blst_p2_in_g2(const POINTonE2& p)
    {
        return (int)POINTonE2_in_G2(p);
    }

    BLST_HD BLST_INLINE int blst_p2_affine_in_g2(const POINTonE2_affine& p)
    {
        POINTonE2 P;

        vec_copy(&P.X, &p.X, 2 * sizeof(P.X));
        vec_select(&P.Z, &p.X, get_BLST_Rx().p, sizeof(P.Z),
            vec_is_zero(&p, sizeof(p)));

        return (int)POINTonE2_in_G2(P);
    }


}

//BLST_EXTERN_C_END
#endif /* BLST_E2_CU_H */