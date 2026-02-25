#ifndef BLST_BACKEND_AVX512_AVX_LAYOUT_H
#define BLST_BACKEND_AVX512_AVX_LAYOUT_H

#include "blst/backend/avx512/avx_ops.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * Shared AVX layout helper counts for radix-48 Fp kernels.
     *
     * Lane notation used below follows avxbls comments:
     * - a lane vector is written as `<H|G|F|E|D|C|B|A>`
     * - each helper applies the same lane rewrite to every vector word in a layout.
     */
    constexpr u64 AVX_FP_NWORDS = 8;
    constexpr u64 AVX_FP_VWORDS = 4;

    /**
     * Per-word lane shuffle: `<H|G|F|E|D|C|B|A> -> <G|G|E|E|C|C|A|A>`.
     * Used to duplicate "component 0" lanes inside an Fp2 pair.
     */
    inline void shuf_00(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VSHUF<0x44>(a[i]);
        }
    }

    /**
     * Per-word lane shuffle: `<H|G|F|E|D|C|B|A> -> <G|H|E|F|C|D|A|B>`.
     * Swaps each adjacent lane-pair.
     */
    inline void shuf_01(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VSHUF<0x4E>(a[i]);
        }
    }

    /**
     * Per-word lane shuffle: `<H|G|F|E|D|C|B|A> -> <H|H|F|F|D|D|B|B>`.
     * Used to duplicate "component 1" lanes inside an Fp2 pair.
     */
    inline void shuf_11(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VSHUF<0xEE>(a[i]);
        }
    }

    /**
     * Zero-masked shuffle: `<H|G|F|E|D|C|B|A> -> <0|H|0|F|0|D|0|B>`.
     * Keeps only the high lane from each adjacent pair.
     */
    inline void shuf_z1(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VZSHUF<0xEE>(0x3333, a[i]);
        }
    }

    /**
     * Per-word blend with mask `0x33`:
     * `<H|G|F|E|D|C|B|A>` with `<P|O|N|M|L|K|J|I>`
     * -> `<H|G|N|M|D|C|J|I>`.
     */
    inline void blend_0x33(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VMBLEND(0x33, a[i], b[i]);
        }
    }

    /**
     * Per-word blend with mask `0x11`:
     * `<H|G|F|E|D|C|B|A>` with `<P|O|N|M|L|K|J|I>`
     * -> `<H|G|F|M|D|C|B|I>`.
     */
    inline void blend_0x11(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VMBLEND(0x11, a[i], b[i]);
        }
    }

    /**
     * Per-word blend with mask `0x55`:
     * `<H|G|F|E|D|C|B|A>` with `<P|O|N|M|L|K|J|I>`
     * -> `<H|O|F|M|D|K|B|I>`.
     */
    inline void blend_0x55(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VMBLEND(0x55, a[i], b[i]);
        }
    }

    /**
     * Per-word blend with mask `0x77`:
     * `<H|G|F|E|D|C|B|A>` with `<P|O|N|M|L|K|J|I>`
     * -> `<H|O|N|M|D|K|J|I>`.
     */
    inline void blend_0x77(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VMBLEND(0x77, a[i], b[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <G|H|H|G|C|D|D|C>`.
     */
    inline void perm_2332(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0xBE>(a[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <H|H|G|G|D|D|C|C>`.
     */
    inline void perm_3322(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0xFA>(a[i]);
        }
    }

    /**
     * Zero-masked permutation: `<H|G|F|E|D|C|B|A> -> <F|E|0|0|B|A|0|0>`.
     */
    inline void perm_10zz(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VZPERM<0x44>(0xCC, a[i]);
        }
    }

    /**
     * Per-word permutation:
     * `<H|G|F|E|D|C|B|A> -> <F|E|F|E|D|C|D|C>`.
     */
    inline void perm_1010(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x44>(a[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <F|E|H|G|B|A|D|C>`.
     */
    inline void perm_1032(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x4E>(a[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <E|F|G|H|A|B|C|D>`.
     */
    inline void perm_0123(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x1B>(a[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <E|H|E|E|A|D|A|A>`.
     */
    inline void perm_0300(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x30>(a[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <F|G|E|E|B|C|A|A>`.
     */
    inline void perm_1200(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x60>(a[i]);
        }
    }

    /**
     * Per-word permutation: `<H|G|F|E|D|C|B|A> -> <G|E|E|E|C|A|A|A>`.
     */
    inline void perm_2000(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x80>(a[i]);
        }
    }

    /**
     * Zero-masked permutation:
     * `<H|G|F|E|D|C|B|A> -> <0|0|H|G|0|0|D|C>`.
     */
    inline void perm_zz32(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_NWORDS; ++i) {
            r[i] = VZPERM<0xEE>(0x33, a[i]);
        }
    }

    /**
     * Half-length (`VWORDS`) version of `perm_1010`.
     * Applies `<H|G|F|E|D|C|B|A> -> <F|E|F|E|D|C|D|C>` on each word.
     */
    inline void perm_1010_hl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_VWORDS; ++i) {
            r[i] = VPERM<0x44>(a[i]);
        }
    }

    /**
     * Half-length (`VWORDS`) version of `perm_1032`.
     */
    inline void perm_1032_hl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_VWORDS; ++i) {
            r[i] = VPERM<0x4E>(a[i]);
        }
    }

    /**
     * Half-length (`VWORDS`) permutation:
     * `<H|G|F|E|D|C|B|A> -> <H|G|H|G|D|C|D|C>`.
     */
    inline void perm_3232_hl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_VWORDS; ++i) {
            r[i] = VPERM<0xEE>(a[i]);
        }
    }

    /**
     * Half-length (`VWORDS`) zero-masked permutation:
     * `<H|G|F|E|D|C|B|A> -> <0|0|H|G|0|0|D|C>`.
     */
    inline void perm_zz32_hl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < AVX_FP_VWORDS; ++i) {
            r[i] = VZPERM<0xEE>(0x33, a[i]);
        }
    }

    /**
     * Half-length (`VWORDS`) variable-lane permutation.
     * `idx` is a per-lane source selector (mod 8) applied to each vector word.
     */
    inline void perm_var_hl(cc512i* r, const cc512i* a, const cc512i idx)
    {
        for (u64 i = 0; i < AVX_FP_VWORDS; ++i) {
            r[i] = VPERMV(idx, a[i]);
        }
    }

    /**
     * Half-length (`VWORDS`) blend with mask `0x0F`:
     * keeps high half from `a` and low half from `b`.
     */
    inline void blend_0x0F_hl(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < AVX_FP_VWORDS; ++i) {
            r[i] = VMBLEND(0x0F, a[i], b[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) version of `shuf_01`.
     */
    inline void shuf_01_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VSHUF<0x4E>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) shuffle:
     * `<H|G|F|E|D|C|B|A> -> <G|G|E|E|C|C|A|A>`.
     */
    inline void shuf_00_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VSHUF<0x44>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) shuffle:
     * `<H|G|F|E|D|C|B|A> -> <H|H|F|F|D|D|B|B>`.
     */
    inline void shuf_11_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VSHUF<0xEE>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) version of `perm_10zz`.
     */
    inline void perm_10zz_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VZPERM<0x44>(0xCC, a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) version of `perm_1032`.
     */
    inline void perm_1032_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x4E>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) version of `perm_zz32`.
     */
    inline void perm_zz32_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VZPERM<0xEE>(0x33, a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) blend with mask `0x33`.
     */
    inline void blend_0x33_dl(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VMBLEND(0x33, a[i], b[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) blend with mask `0x55`.
     */
    inline void blend_0x55_dl(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VMBLEND(0x55, a[i], b[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) permutation:
     * `<H|G|F|E|D|C|B|A> -> <F|E|F|E|D|C|D|C>`.
     */
    inline void perm_1010_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x44>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) permutation:
     * `<H|G|F|E|D|C|B|A> -> <F|E|H|G|B|A|D|C>`.
     */
    inline void perm_1100_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0x50>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) permutation:
     * `<H|G|F|E|D|C|B|A> -> <G|H|G|H|C|D|C|D>`.
     */
    inline void perm_2323_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0xBB>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) permutation:
     * `<H|G|F|E|D|C|B|A> -> <H|E|E|E|D|A|A|A>`.
     */
    inline void perm_3000_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0xC0>(a[i]);
        }
    }

    /**
     * Double-length (`2*NWORDS`) permutation:
     * `<H|G|F|E|D|C|B|A> -> <H|G|H|G|D|C|D|C>`.
     */
    inline void perm_3232_dl(cc512i* r, const cc512i* a)
    {
        for (u64 i = 0; i < 2 * AVX_FP_NWORDS; ++i) {
            r[i] = VPERM<0xEE>(a[i]);
        }
    }

    /**
     * One-and-half-length (`3*VWORDS`) variable permutation.
     * Used by widened `4x2w`/`2x2x2w` kernels.
     */
    inline void perm_var_vl(cc512i* r, const cc512i* a, const cc512i idx)
    {
        for (u64 i = 0; i < 3 * AVX_FP_VWORDS; ++i) {
            r[i] = VPERMV(idx, a[i]);
        }
    }

    /**
     * One-and-half-length (`3*VWORDS`) blend with mask `0x0F`.
     * Used by widened `1x2x2x2w`/`2x2x2x1w` kernels.
     */
    inline void blend_0x0F_vl(cc512i* r, const cc512i* a, const cc512i* b)
    {
        for (u64 i = 0; i < 3 * AVX_FP_VWORDS; ++i) {
            r[i] = VMBLEND(0x0F, a[i], b[i]);
        }
    }
}

#endif

#endif // BLST_BACKEND_AVX512_AVX_LAYOUT_H
