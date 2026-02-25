#ifndef BLST_FIELDS_FP_AVX_H
#define BLST_FIELDS_FP_AVX_H

#include "blst/common/defines.h"
#include "blst/backend/avx512/avx_ops.h"
#include <cstddef>
#include <type_traits>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 Fp helpers ported from avxbls.
     *
     * Representation notes:
     * - Field arithmetic is over BLS12-381 Fp in Montgomery form.
     * - Limb radix is 2^48. Reduced values use 8 radix-48 limbs.
     * - Most routines operate on packed vectors (multiple field elements in parallel).
     * - Layout-specific helpers intentionally use avxbls naming to keep ported kernels readable.
     */

    /** @name AVX radix-48 parameters (BLS12-381 Fp) */
    ///@{
    /// @brief Bit-width of one internal limb in AVX kernels.
    constexpr u64 BRADIX = 48;
    /// @brief Number of radix-48 limbs in one single-length Fp value.
    constexpr u64 NWORDS = 8;
    /// @brief Number of vector words in the `4x2w` packing.
    constexpr u64 VWORDS = 4;
    /// @brief Number of 64-bit limbs in canonical scalar Fp storage.
    constexpr u64 SWORDS = 6;
    /// @brief Mask for one radix-48 limb.
    constexpr u64 BMASK = 0xFFFFFFFFFFFFULL;
    /// @brief Mask for one full 64-bit word.
    constexpr u64 FMASK = 0xFFFFFFFFFFFFFFFFULL;
    /// @brief Left-shift used to align MADD52 high partials into radix-48 limbs.
    constexpr u64 BALIGN = 4;
    /// @brief Montgomery reduction constant for radix-48 kernels (`-p^{-1} mod 2^48`).
    constexpr u64 MONT_W_R48 = 0xFFFCFFFCFFFDULL;
    ///@}

    /// @brief Modulus in radix-48 limbs (little-endian limb order).
    extern const u64 P48[NWORDS];

    template <std::size_t N>
    struct AvxPack
    {
        cc512i lanes[N];

        BLST_FORCEINLINE cc512i& operator[](std::size_t i)
        {
            BLST_ASSERT(i < N);
            return lanes[i];
        }

        BLST_FORCEINLINE const cc512i& operator[](std::size_t i) const
        {
            BLST_ASSERT(i < N);
            return lanes[i];
        }

        BLST_FORCEINLINE cc512i* data()
        {
            return lanes;
        }

        BLST_FORCEINLINE const cc512i* data() const
        {
            return lanes;
        }

        BLST_FORCEINLINE operator cc512i*()
        {
            return lanes;
        }

        BLST_FORCEINLINE operator const cc512i*() const
        {
            return lanes;
        }
    };

    static_assert(std::is_trivially_default_constructible_v<AvxPack<1>>);
    static_assert(std::is_trivially_copyable_v<AvxPack<1>>);

#if defined(NDEBUG)
#define BLST_ASSERT_NO_OVERLAP_AVX(a, b) ((void)0)
#else
    BLST_FORCEINLINE bool avx_ranges_overlap(const void* a, std::size_t a_bytes, const void* b, std::size_t b_bytes)
    {
        const auto* ap = static_cast<const u8*>(a);
        const auto* bp = static_cast<const u8*>(b);
        return (ap < (bp + b_bytes)) && (bp < (ap + a_bytes));
    }

    template <std::size_t NA, std::size_t NB>
    BLST_FORCEINLINE void avx_assert_no_overlap(const AvxPack<NA>& a, const AvxPack<NB>& b)
    {
        BLST_ASSERT(!avx_ranges_overlap(a.data(), sizeof(a), b.data(), sizeof(b)));
    }

#define BLST_ASSERT_NO_OVERLAP_AVX(a, b) ::blst::avx_assert_no_overlap((a), (b))
#endif

    /** @name Packed AVX Fp layouts */
    ///@{
    /// @brief Eight parallel reduced Fp values in "8x1w" layout.
    /// @details
    /// Logical mapping:
    /// - `value[ch].limb[i] == a[i].lane[ch]`
    /// - `i in [0,7]`, `ch in [0,7]`
    /// - each limb is expected in `[0, 2^48)`.
    using FpPack8x1 = AvxPack<NWORDS>;
    using fp_8x1w = FpPack8x1;

    /// @brief Widened container for `fp_8x1w` products and linear combinations.
    /// @details
    /// Stores 16 radix-48 limbs per channel (unreduced double-length value).
    /// Logical mapping is the same as `fp_8x1w` but over indices `i in [0,15]`.
    using UFpPack8x1 = AvxPack<2 * NWORDS>;
    using fpx2_8x1w = UFpPack8x1;

    /// @brief Four parallel reduced Fp values in "4x2w" lane-pair layout.
    /// @details
    /// For element `e in [0,3]` and limb `j in [0,7]`:
    /// - `word = j % 4`
    /// - `lane = 2 * e + (j / 4)` (`j/4` is 0 for low half, 1 for high half)
    /// - `value[e].limb[j] == a[word].lane[lane]`
    /// Each limb is expected in `[0, 2^48)`.
    using FpPack4x2 = AvxPack<VWORDS>;
    using fp_4x2w = FpPack4x2;

    /// @brief Widened container paired with `fp_4x2w` kernels.
    /// @details
    /// Contains 12 vector words (`3 * VWORDS`) as used by avxbls 4x2w double-length kernels.
    /// This is not a standalone canonical scalar layout; interpret via matching kernels.
    using UFpPack4x2 = AvxPack<3 * VWORDS>;
    using fpx2_4x2w = UFpPack4x2;
    ///@}

    static_assert(std::is_trivially_default_constructible_v<FpPack8x1>);
    static_assert(std::is_trivially_copyable_v<FpPack8x1>);
    static_assert(std::is_trivially_default_constructible_v<UFpPack8x1>);
    static_assert(std::is_trivially_copyable_v<UFpPack8x1>);
    static_assert(std::is_trivially_default_constructible_v<FpPack4x2>);
    static_assert(std::is_trivially_copyable_v<FpPack4x2>);
    static_assert(std::is_trivially_default_constructible_v<UFpPack4x2>);
    static_assert(std::is_trivially_copyable_v<UFpPack4x2>);

    /** @name Conversion and layout helpers */
    ///@{
    /// @brief Convert one scalar from 6x64-bit limbs to 8x48-bit limbs.
    /// @pre `a` points to at least 6 limbs; `r` points to at least 8 limbs.
    /// @post `r[0..7]` stores the same integer split into radix-48 limbs (`r[i] < 2^48`).
    void conv_64to48_fp(u64* r, const u64* a);

    /// @brief Lane-wise `conv_64to48_fp` for 8 packed channels.
    /// @pre `a` has 6 vectors of 64-bit limbs; `r` has 8 vectors for radix-48 limbs.
    /// @post Each lane `ch` is converted exactly as `conv_64to48_fp`.
    void conv_64to48_fp_8x1w(fp_8x1w& r, const cc512i* BLST_restrict a);

    /// @brief Lane-pair `conv_64to48_fp` for 4 packed channels in 4x2w layout.
    /// @pre `a` has 3 vectors (64-bit packing), `r` has 4 vectors (radix-48 packing).
    /// @post Each packed channel is converted into 8 radix-48 limbs.
    void conv_64to48_fp_4x2w(fp_4x2w& r, const cc512i* BLST_restrict a);

    /// @brief Generic multi-precision conversion from base 2^64 to base 2^48.
    /// @param rlen Number of output radix-48 limbs.
    /// @param alen Number of input 64-bit limbs.
    /// @post High unused output limbs are zero-filled.
    void conv_64to48_mpi(u64* r, const u64* a, int rlen, int alen);


    /// @brief Convert one scalar from 8x48-bit limbs to 6x64-bit limbs.
    /// @pre `a` points to 8 radix-48 limbs.
    /// @post `r` reconstructs the same integer in 64-bit limb form.
    void conv_48to64_fp(u64* r, const u64* a);

    /// @brief Lane-wise `conv_48to64_fp` for 8 packed channels.
    void conv_48to64_fp_8x1w(cc512i* BLST_restrict r, const fp_8x1w& a);

    /// @brief Lane-pair `conv_48to64_fp` for 4 packed channels in 4x2w layout.
    void conv_48to64_fp_4x2w(cc512i* BLST_restrict r, const fp_4x2w& a);

    /// @brief Generic multi-precision conversion from base 2^48 to base 2^64.
    /// @post High unused output limbs are zero-filled.
    void conv_48to64_mpi(u64* r, const u64* a, int rlen, int alen);


    /// @brief Carry-propagate scalar radix-48 limbs in place.
    /// @pre `a` has at least `NWORDS` limbs.
    /// @post `a[0..NWORDS-2] < 2^48`; final carry remains in `a[NWORDS-1]`.
    void carryp_mpi48(u64* a);

    /// @brief Carry-propagate packed `fp_4x2w` limbs in place.
    /// @post Low three limb words are masked to 48 bits per lane.
    void carryp_fp_4x2w(fp_4x2w& a);

    /// @brief Extract one `8x1w` lane into scalar radix-48 limbs.
    /// @param ch Lane index in `[0,7]`.
    /// @post `r[i] = a[i].lane[ch]` for `i in [0,7]`.
    void get_channel_8x1w(u64* BLST_restrict r, const fp_8x1w& a, int ch);

    /// @brief Extract one `4x2w` packed channel (two adjacent lanes).
    /// @param ch Base lane index (typically even), valid range `[0,6]`.
    /// @post `r[0..3]` come from lane `ch`, `r[4..7]` from lane `ch+1`.
    void get_channel_4x2w(u64* BLST_restrict r, const fp_4x2w& a, int ch);
    ///@}



    /** @name Reduced Fp helpers */
    ///@{
    /// @brief Lane-wise modular add for reduced `fp_8x1w`.
    /// @pre Each lane of `a` and `b` represents a reduced Montgomery Fp value.
    /// @post `r[ch] = a[ch] + b[ch] (mod p)` and each output limb is masked to 48 bits.
    void add_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);

    /// @brief Lane-wise modular subtract for reduced `fp_8x1w`.
    /// @post `r[ch] = a[ch] - b[ch] (mod p)`, reduced to radix-48 limbs.
    void sub_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);

    /// @brief Mixed add/sub helper for reduced `fp_8x1w` with alternating lane schedule.
    /// @details With lane order `<H|G|F|E|D|C|B|A>` and `<P|O|N|M|L|K|J|I>`:
    ///          computes `<H+P, G-O, F+N, E-M, D+L, C-K, B+J, A-I>` modulo `p`.
    ///          Equivalently by lane index: even lanes subtract, odd lanes add.
    /// @post Result is reduced modulo `p` per lane.
    void asx4_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);

    /// @brief Mixed add/sub helper for reduced `fp_8x1w` (`asax6` schedule).
    /// @details Computes `<H+P, G-O, F+N, E+M, D+L, C+K, B+J, A+I>` modulo `p`.
    ///          Lane rule: add on lanes 0..5 and 7, subtract on lane 6.
    void asax6_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);

    /// @brief Mixed add/sub helper for reduced `fp_8x1w` (`assx6` schedule).
    /// @details Computes `<H+P, G-O, F-N, E-M, D-L, C-K, B-J, A-I>` modulo `p`.
    ///          Lane rule: subtract on lanes 0..6, add on lane 7.
    void assx6_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);


    /// @brief Lane-wise modular add for reduced `fp_4x2w`.
    /// @post Each packed channel equals `a + b (mod p)` in 4x2w layout.
    void add_fp_4x2w(fp_4x2w& r, const fp_4x2w& a, const fp_4x2w& b);

    /// @brief Lane-wise modular subtract for reduced `fp_4x2w`.
    /// @post Each packed channel equals `a - b (mod p)` in 4x2w layout.
    void sub_fp_4x2w(fp_4x2w& r, const fp_4x2w& a, const fp_4x2w& b);

    /// @brief Mixed add/sub helper for reduced `fp_4x2w` (`asx2` schedule).
    /// @details Using avxbls symbols:
    ///          `<D'+H'-p, D+H-p, C'-G', C-G, B'+F'-p', B+F-p, A'-E', A-E>` then final conditional add.
    /// @post Reduced modulo `p` in 4x2w layout.
    void asx2_fp_4x2w(fp_4x2w& r, const fp_4x2w& a, const fp_4x2w& b);
    ///@}


    /** @name Widened Fp helper combinations */
    ///@{
    /// @brief Lane-wise add on widened `fpx2_8x1w` values.
    /// @post `r = a + b` with carry propagation and conditional correction by `p * 2^384`.
    void add_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);

    /// @brief Lane-wise subtract on widened `fpx2_8x1w` values.
    /// @post `r = a - b` with carry propagation and conditional correction by `p * 2^384`.
    void sub_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);

    /// @brief Widened mixed helper (`aaasx2`) on `fpx2_8x1w`.
    /// @details Computes `<H+P, G+O, F+N, E-M, D+L, C+K, B+J, A-I>` on low half,
    ///          with matching high-half correction (`-p` on add lanes), then normalizes.
    void aaasx2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);

    /// @brief Widened mixed helper (`aassx2`) on `fpx2_8x1w`.
    /// @details Computes `<H+P, G+O, F-N, E-M, D+L, C+K, B-J, A-I>` with high-half correction and normalization.
    void aassx2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);

    /// @brief Widened mixed helper (`asx4`) on `fpx2_8x1w`.
    /// @details Computes `<H+P, G-O, F+N, E-M, D+L, C-K, B+J, A-I>` with high-half correction and normalization.
    void asx4_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);

    /// @brief Widened mixed helper (`ssaax2`) on `fpx2_8x1w`.
    /// @details Computes `<H-P, G-O, F+N, E+M, D-L, C-K, B+J, A+I>` with high-half correction and normalization.
    void ssaax2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);

    /// @brief Widened mixed helper (`ssasx2`) on `fpx2_8x1w`.
    /// @details Computes `<H-P, G-O, F+N, E-M, D-L, C-K, B+J, A-I>` with high-half correction and normalization.
    void ssasx2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b);


    /// @brief Lane-wise add on widened `fpx2_4x2w` values (4x2w schedule).
    /// @post `r = a + b` with carry propagation and conditional correction by packed `p`.
    void add_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b);

    /// @brief Lane-wise subtract on widened `fpx2_4x2w` values (4x2w schedule).
    /// @post `r = a - b` with carry propagation and conditional correction by packed `p`.
    void sub_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b);

    /// @brief Widened mixed helper (`asx2`) on `fpx2_4x2w`.
    /// @details Computes
    /// `<D'+H'-p, D+H-p, C'-G', C-G, B'+F'-p', B+F-p, A'-E', A-E>` then normalizes.
    void asx2_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b);

    /// @brief Widened mixed helper (`aass`) on `fpx2_4x2w`.
    /// @details Computes
    /// `<D'+H'-p, D+H-p, C'+G'-p', C+G-p, B'-F', B-F, A'-E', A-E>` then normalizes.
    void aass_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b);
    ///@}



    /** @name Fp multiply and reduction kernels */
    ///@{
    /// @brief Canonical widened multiply for `fp_8x1w`.
    /// @details This codebase intentionally uses avxbls v2 only.
    /// @post Writes unreduced double-length product to `fpx2_8x1w`.
    void mul_fpx2_8x1w(fpx2_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);

    /// @brief AVXBLS v2 widened multiply kernel for `fp_8x1w`.
    /// @post Writes unreduced widened product into `fpx2_8x1w` (same layout as `mul_fpx2_8x1w`).
    void mul_fpx2_8x1w_v2(fpx2_8x1w& r, const fp_8x1w& a, const fp_8x1w& b);

    /// @brief Canonical widened multiply for `fp_4x2w`.
    /// @details Current implementation forwards to `mul_fpx2_4x2w_v1`.
    /// @post Writes unreduced widened product into `fpx2_4x2w`.
    void mul_fpx2_4x2w(fpx2_4x2w& r, const fp_4x2w& a, const fp_4x2w& b);

    /// @brief AVXBLS widened multiply kernel for `fp_4x2w`.
    /// @post Writes unreduced widened product into `fpx2_4x2w`.
    void mul_fpx2_4x2w_v1(fpx2_4x2w& r, const fp_4x2w& a, const fp_4x2w& b);


    /// @brief Montgomery reduction from widened `8x1w` to reduced `fp_8x1w`.
    /// @pre `a` is a valid double-length radix-48 value (typically from `mul_fpx2_8x1w*`).
    /// @post `r` is reduced modulo `p` in Montgomery form, radix-48 limbs.
    void redc_fpx2_8x1w(fp_8x1w& r, const fpx2_8x1w& a);

    /// @brief Montgomery reduction from widened `4x2w` to reduced `fp_4x2w`.
    /// @pre `a` is a valid double-length 4x2w value (typically from `mul_fpx2_4x2w*`).
    /// @post `r` is reduced modulo `p` in Montgomery form, 4x2w layout.
    void redc_fpx2_4x2w(fp_4x2w& r, const fpx2_4x2w& a);
    ///@}
}

#endif

#endif // BLST_FIELDS_FP_AVX_H
