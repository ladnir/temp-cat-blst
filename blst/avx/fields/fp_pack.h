#ifndef BLST_FIELDS_FP_AVX_H
#define BLST_FIELDS_FP_AVX_H

#include "blst/common/defines.h"
#include "blst/backend/avx512/avx_ops.h"
#include "blst/fields/Fp.h"
#include <cstddef>
#include <type_traits>
#include <array>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#if !defined(BLST_AVX_NOINLINE)
# if defined(_MSC_VER)
#  define BLST_AVX_NOINLINE __declspec(noinline)
# elif defined(__GNUC__) || defined(__clang__)
#  define BLST_AVX_NOINLINE __attribute__((noinline))
# else
#  define BLST_AVX_NOINLINE
# endif
#endif

namespace blst::avx
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

    inline constexpr const std::array< u64, NWORDS> P48 = {
        0xFFFFFFFFAAABULL, 0xB153FFFFB9FEULL, 0xF6241EABFFFEULL, 0x6730D2A0F6B0ULL,
        0x4B84F38512BFULL, 0x434BACD76477ULL, 0xE69A4B1BA7B6ULL, 0x1A0111EA397FULL,
    };


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

#define BLST_ASSERT_NO_OVERLAP_AVX(a, b) ::blst::avx::avx_assert_no_overlap((a), (b))
#endif

    struct LFpPack8x1;
    struct LUFpPack8x1;
    struct UFpPack8x1;

    /// @brief Compile-time renormalization policy for lazy widened `LUFpPack8x1`.
    enum class LUFpRenormClass
    {
        Tight,
        Medium,
        Loose,
    };

    /// @brief Compile-time bound class for fast `LUFpPack8x1 -> UFpPack8x1` normalization.
    enum class LUFpNormalizeBound
    {
        /// @brief After signed carry, the high half satisfies `-p <= H < 2p`.
        HiNeg1Pos2,
        /// @brief Intended for `Fp6 mul` final-repair inputs derived from the
        ///        `LNeg4_Pos4_HNeg4_Pos4` worksheet class.
        /// @details The specialization is deliberately narrow and benchmarked
        ///          only for that exact use case.
        HiNeg4Pos4,
    };

    /** @name Packed AVX Fp layouts */
    ///@{
    /// @brief Eight parallel reduced Fp values in "8x1w" layout.
    /// @details
    /// Logical mapping:
    /// - `value[ch].limb[i] == a[i].lane[ch]`
    /// - `i in [0,7]`, `ch in [0,7]`
    /// - each limb is expected in `[0, 2^48)`.
    struct FpPack8x1 : public AvxPack<NWORDS>
    {

        /// @brief Lane-wise `conv_64to48_fp` for 8 packed channels.
        /// @pre `a` has 6 vectors of 64-bit limbs; `r` has 8 vectors for radix-48 limbs.
        /// @post Each lane `ch` is converted exactly as `conv_64to48_fp`.
        static void conv_64to48(FpPack8x1& r, const cc512i* BLST_restrict a);

        /// @brief Lane-wise `conv_48to64_fp` for 8 packed channels.
        static void conv_48to64(cc512i* BLST_restrict r, const FpPack8x1& a);

        /// @brief Pack 8 scalar `Fp` values (AoS) into one `FpPack8x1`.
        /// @pre `in8` points to at least 8 reduced Montgomery `Fp` values.
        /// @post lane `i` of `r` equals `in8[i]` in radix-48 representation.
        static void from_aos(FpPack8x1& r, const Fp* in8);

        /// @brief Unpack one `FpPack8x1` to 8 scalar `Fp` values (AoS).
        /// @pre `out8` points to at least 8 `Fp` values.
        /// @post `out8[i]` equals lane `i` of `a`.
        static void to_aos(Fp* out8, const FpPack8x1& a);

        /// @brief Partial AoS pack for tail handling.
        /// @param n_active Number of active lanes in `[0,8]`.
        /// @param fill Value used for inactive lanes.
        /// @post lanes `[0,n_active)` are loaded from `in`, lanes `[n_active,8)` are `fill`.
        static void from_aos_partial(FpPack8x1& r, const Fp* in, std::size_t n_active, const Fp& fill);

        /// @brief Partial AoS unpack for tail handling.
        /// @param n_active Number of active lanes in `[0,8]`.
        /// @post only `out[0..n_active)` are written.
        static void to_aos_partial(Fp* out, std::size_t n_active, const FpPack8x1& a);

        /// @brief Lane-wise select between two packed values.
        /// @details Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b, u8 lane_mask);



        /** @name Reduced Fp helpers */
        ///@{
        /// @brief Lane-wise modular add for reduced `FpPack8x1`.
        /// @pre Each lane of `a` and `b` represents a reduced Montgomery Fp value.
        /// @post `r[ch] = a[ch] + b[ch] (mod p)` and each output limb is masked to 48 bits.
        static void add(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);

        /// @brief Lane-wise modular subtract for reduced `FpPack8x1`.
        /// @post `r[ch] = a[ch] - b[ch] (mod p)`, reduced to radix-48 limbs.
        static void sub(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);

        /// @brief Mixed add/sub helper for reduced `FpPack8x1` with alternating lane schedule.
        /// @details With lane order `<H|G|F|E|D|C|B|A>` and `<P|O|N|M|L|K|J|I>`:
        ///          computes `<H+P, G-O, F+N, E-M, D+L, C-K, B+J, A-I>` modulo `p`.
        ///          Equivalently by lane index: even lanes subtract, odd lanes add.
        /// @post Result is reduced modulo `p` per lane.
        static void asx4(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);

        /// @brief Mixed add/sub helper for reduced `FpPack8x1` (`asax6` schedule).
        /// @details Computes `<H+P, G-O, F+N, E+M, D+L, C+K, B+J, A+I>` modulo `p`.
        ///          Lane rule: add on lanes 0..5 and 7, subtract on lane 6.
        static void asax6(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);

        /// @brief Mixed add/sub helper for reduced `FpPack8x1` (`assx6` schedule).
        /// @details Computes `<H+P, G-O, F-N, E-M, D-L, C-K, B-J, A-I>` modulo `p`.
        ///          Lane rule: subtract on lanes 0..6, add on lane 7.
        static void assx6(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);

        /// @brief Extract one `8x1w` lane into scalar radix-48 limbs.
        /// @param ch Lane index in `[0,7]`.
        /// @post `r[i] = a[i].lane[ch]` for `i in [0,7]`.
        static void get_channel(u64* BLST_restrict r, const FpPack8x1& a, int ch);


        /// @brief Canonical widened multiply for `FpPack8x1`.
        /// @post Writes unreduced double-length product to `UFpPack8x1`.
        static void mul(UFpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);


    };

    /// @brief Eight parallel lazy-reduced Fp values in "8x1w" layout.
    /// @details
    /// Same physical layout as `FpPack8x1`, but semantic range is intentionally
    /// non-canonical (`[0, c_l * p)` for a small constant `c_l`).
    struct LFpPack8x1 : public AvxPack<NWORDS>
    {
    };

    /// @brief Lazy widened companion for `UFpPack8x1` in "8x1w" layout.
    /// @details
    /// Same physical 16-limb radix-48 storage as `UFpPack8x1`, but semantic range
    /// is intentionally relaxed and signed for local linear-combination regions.
    ///
    /// Representation model:
    /// - one widened value per AVX lane
    /// - value uses radix `R = 2^48`
    /// - limbs may be negative and may temporarily exceed 48 bits
    ///
    /// Contract:
    /// - `LUFpPack8x1` is not accepted by existing `UFpPack8x1::redc_core`
    /// - `renorm<B>` keeps the value in the lazy widened domain while shrinking
    ///   limb growth
    /// - `normalize(UFpPack8x1&, const LUFpPack8x1&)` restores the strict widened
    ///   `UFpPack8x1` contract at an explicit boundary
    struct LUFpPack8x1 : public AvxPack<2 * NWORDS>
    {
        /// @brief Lane-wise lazy widened add.
        /// @post `r = a + b` in the lazy widened signed domain, without restoring
        ///       the strict `UFpPack8x1` bounded contract.
        static void add(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b);
        static void add(LUFpPack8x1& r, const UFpPack8x1& a, const LUFpPack8x1& b);
        static void add(LUFpPack8x1& r, const LUFpPack8x1& a, const UFpPack8x1& b);
        static void add(LUFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Lane-wise lazy widened subtract.
        /// @post `r = a - b` in the lazy widened signed domain, without restoring
        ///       the strict `UFpPack8x1` bounded contract.
        static void sub(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b);
        static void sub(LUFpPack8x1& r, const UFpPack8x1& a, const LUFpPack8x1& b);
        static void sub(LUFpPack8x1& r, const LUFpPack8x1& a, const UFpPack8x1& b);
        static void sub(LUFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Compile-time renormalization for the lazy widened domain.
        /// @details `B` selects the expected drift class known by the caller.
        ///          The represented integer value is unchanged.
        template <LUFpRenormClass B>
        static void renorm(LUFpPack8x1& r, const LUFpPack8x1& a);

        /// @brief Normalize lazy widened storage back into strict widened `UFpPack8x1`.
        /// @details This is an explicit boundary operation. It performs any
        ///          required carry propagation and high-half correction so the
        ///          output again satisfies the current `UFpPack8x1` contract.
        static void normalize(UFpPack8x1& r, const LUFpPack8x1& a);

        /// @brief Fast bounded normalize from lazy widened to strict widened storage.
        /// @details The caller is responsible for maintaining the documented
        ///          compile-time bound class for the high half after signed carry.
        template <LUFpNormalizeBound B>
        static void normalize(UFpPack8x1& r, const LUFpPack8x1& a);

        /// @brief Lane-wise select between two lazy widened packed values.
        /// @details Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b, u8 lane_mask);

    private:
        static void normalize_hi_neg1_pos2(UFpPack8x1& r, const LUFpPack8x1& a);
        static void normalize_hi_neg4_pos4(UFpPack8x1& r, const LUFpPack8x1& a);
    };

    /// @brief Redundant widened companion for `UFpPack8x1` in "8x1w" layout.
    /// @details
    /// Same physical 16-limb radix-48 storage as `UFpPack8x1`, but intended for
    /// fully carry-lazy local combine regions.
    ///
    /// Representation model:
    /// - one widened value per AVX lane
    /// - value uses radix `R = 2^48`
    /// - limbs may be negative and may temporarily exceed 48 bits
    /// - no inter-limb carry propagation is implied by `add/sub`
    /// @brief Widened container for `FpPack8x1` products and linear combinations.
    /// @details
    /// Stores 16 radix-48 limbs per channel (unreduced double-length value).
    /// Logical mapping is the same as `FpPack8x1` but over indices `i in [0,15]`.
    struct UFpPack8x1 : public  AvxPack<2 * NWORDS>
    {

        /// @brief Lane-wise add on widened `UFpPack8x1` values.
        /// @post `r = a + b` with carry propagation and conditional correction by `p * 2^384`.
        static void add(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Lane-wise subtract on widened `UFpPack8x1` values.
        /// @post `r = a - b` with carry propagation and conditional correction by `p * 2^384`.
        static void sub(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Widened mixed helper (`aaasx2`) on `UFpPack8x1`.
        /// @details Computes `<H+P, G+O, F+N, E-M, D+L, C+K, B+J, A-I>` on low half,
        ///          with matching high-half correction (`-p` on add lanes), then normalizes.
        static void aaasx2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Widened mixed helper (`aassx2`) on `UFpPack8x1`.
        /// @details Computes `<H+P, G+O, F-N, E-M, D+L, C+K, B-J, A-I>` with high-half correction and normalization.
        static void aassx2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Widened mixed helper (`asx4`) on `UFpPack8x1`.
        /// @details Computes `<H+P, G-O, F+N, E-M, D+L, C-K, B+J, A-I>` with high-half correction and normalization.
        static void asx4(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Widened mixed helper (`ssaax2`) on `UFpPack8x1`.
        /// @details Computes `<H-P, G-O, F+N, E+M, D-L, C-K, B+J, A+I>` with high-half correction and normalization.
        static void ssaax2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Widened mixed helper (`ssasx2`) on `UFpPack8x1`.
        /// @details Computes `<H-P, G-O, F+N, E-M, D-L, C-K, B+J, A-I>` with high-half correction and normalization.
        static void ssasx2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b);

        /// @brief Montgomery reduction core from widened `8x1w` to lazy-reduced `LFpPack8x1`.
        /// @pre `a` is a valid double-length radix-48 value (typically from `FpPack8x1::mul*`).
        /// @post `r` has carry-propagated radix-48 limbs without final canonical subtraction.
        static void redc_core(LFpPack8x1& r, const UFpPack8x1& a);

        /// @brief Canonicalize a lazy-reduced `LFpPack8x1` into reduced `FpPack8x1`.
        /// @pre `a` is a valid lazy-reduced radix-48 value.
        /// @post `r` is reduced modulo `p` in Montgomery form, radix-48 limbs.
        static void normalize(FpPack8x1& r, const LFpPack8x1& a);

        /// @brief Montgomery reduction from widened `8x1w` to reduced `FpPack8x1`.
        /// @pre `a` is a valid double-length radix-48 value (typically from `FpPack8x1::mul*`).
        /// @post `r` is reduced modulo `p` in Montgomery form, radix-48 limbs.
        static void redc(FpPack8x1& r, const UFpPack8x1& a);

        /// @brief Lane-wise select between two widened packed values.
        /// @details Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b, u8 lane_mask);
    };

    struct UFpPack4x2;

    /// @brief Four parallel reduced Fp values in "4x2w" lane-pair layout.
    /// @details
    /// For element `e in [0,3]` and limb `j in [0,7]`:
    /// - `word = j % 4`
    /// - `lane = 2 * e + (j / 4)` (`j/4` is 0 for low half, 1 for high half)
    /// - `value[e].limb[j] == a[word].lane[lane]`
    /// Each limb is expected in `[0, 2^48)`.
    struct FpPack4x2 : public AvxPack<VWORDS> {

        /// @brief Lane-pair `conv_64to48_fp` for 4 packed channels in 4x2w layout.
        /// @pre `a` has 3 vectors (64-bit packing), `r` has 4 vectors (radix-48 packing).
        /// @post Each packed channel is converted into 8 radix-48 limbs.
        static void conv_64to48(FpPack4x2& r, const cc512i* BLST_restrict a);


        /// @brief Lane-pair `conv_48to64_fp` for 4 packed channels in 4x2w layout.
        static void conv_48to64(cc512i* BLST_restrict r, const FpPack4x2& a);


        /// @brief Carry-propagate packed `FpPack4x2` limbs in place.
        /// @post Low three limb words are masked to 48 bits per lane.
        static void carryp(FpPack4x2& a);

        /// @brief Lane-wise modular add for reduced `FpPack4x2`.
        /// @post Each packed channel equals `a + b (mod p)` in 4x2w layout.
        static void add(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);

        /// @brief Lane-wise modular subtract for reduced `FpPack4x2`.
        /// @post Each packed channel equals `a - b (mod p)` in 4x2w layout.
        static void sub(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);

        /// @brief Mixed add/sub helper for reduced `FpPack4x2` (`asx2` schedule).
        /// @details Using avxbls symbols:
        ///          `<D'+H'-p, D+H-p, C'-G', C-G, B'+F'-p', B+F-p, A'-E', A-E>` then final conditional add.
        /// @post Reduced modulo `p` in 4x2w layout.
        static void asx2(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);


        /// @brief Canonical widened multiply for `FpPack4x2`.
        /// @details Current implementation forwards to `mul_v1`.
        /// @post Writes unreduced widened product into `UFpPack4x2`.
        static void mul(UFpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);

        /// @brief AVXBLS widened multiply kernel for `FpPack4x2`.
        /// @post Writes unreduced widened product into `UFpPack4x2`.
        static void mul_v1(UFpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);



        /// @brief Extract one `4x2w` packed channel (two adjacent lanes).
        /// @param ch Base lane index (typically even), valid range `[0,6]`.
        /// @post `r[0..3]` come from lane `ch`, `r[4..7]` from lane `ch+1`.
        static void get_channel(u64* BLST_restrict r, const FpPack4x2& a, int ch);
    };
    //using FpPack4x2 = FpPack4x2;

    /// @brief Widened container paired with `FpPack4x2` kernels.
    /// @details
    /// Contains 12 vector words (`3 * VWORDS`) as used by avxbls 4x2w double-length kernels.
    /// This is not a standalone canonical scalar layout; interpret via matching kernels.
    struct UFpPack4x2 : AvxPack<3 * VWORDS> {
    


        /// @brief Lane-wise add on widened `UFpPack4x2` values (4x2w schedule).
        /// @post `r = a + b` with carry propagation and conditional correction by packed `p`.
        static void add(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b);

        /// @brief Lane-wise subtract on widened `UFpPack4x2` values (4x2w schedule).
        /// @post `r = a - b` with carry propagation and conditional correction by packed `p`.
        static void sub(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b);

        /// @brief Widened mixed helper (`asx2`) on `UFpPack4x2`.
        /// @details Computes
        /// `<D'+H'-p, D+H-p, C'-G', C-G, B'+F'-p', B+F-p, A'-E', A-E>` then normalizes.
        static void asx2(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b);

        /// @brief Widened mixed helper (`aass`) on `UFpPack4x2`.
        /// @details Computes
        /// `<D'+H'-p, D+H-p, C'+G'-p', C+G-p, B'-F', B-F, A'-E', A-E>` then normalizes.
        static void aass(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b);

        /// @brief Montgomery reduction from widened `4x2w` to reduced `FpPack4x2`.
        /// @pre `a` is a valid double-length 4x2w value (typically from `FpPack4x2::mul*`).
        /// @post `r` is reduced modulo `p` in Montgomery form, 4x2w layout.
        static void redc(FpPack4x2& r, const UFpPack4x2& a);
    
    };

    template <LUFpRenormClass B>
    BLST_FORCEINLINE void LUFpPack8x1::renorm(LUFpPack8x1& r, const LUFpPack8x1& a)
    {
        BLST_ASSERT_NO_OVERLAP_AVX(r, a);

        cc512i limbs[2 * NWORDS];
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            limbs[i] = a[i];
        }

        for (std::size_t i = 0; i < (2 * NWORDS) - 1; ++i) {
            const cc512i carry = VSRA(limbs[i], BRADIX);
            limbs[i] = VSUB(limbs[i], VSHL(carry, BRADIX));
            limbs[i + 1] = VADD(limbs[i + 1], carry);
        }

        if constexpr (B == LUFpRenormClass::Tight) {
            // One signed carry pass is sufficient for the initial conservative implementation.
        } else if constexpr (B == LUFpRenormClass::Medium) {
            // Reserved for future tighter medium-drift repair once bounds are tuned.
        } else {
            static_assert(B == LUFpRenormClass::Loose);
            // Reserved for future loose-drift repair once bounds are tuned.
        }

        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            r[i] = limbs[i];
        }
    }

    template <LUFpNormalizeBound B>
    BLST_FORCEINLINE void LUFpPack8x1::normalize(UFpPack8x1& r, const LUFpPack8x1& a)
    {
        if constexpr (B == LUFpNormalizeBound::HiNeg1Pos2) {
            normalize_hi_neg1_pos2(r, a);
        } else if constexpr (B == LUFpNormalizeBound::HiNeg4Pos4) {
            normalize_hi_neg4_pos4(r, a);
        } else {
            static_assert(
                B == LUFpNormalizeBound::HiNeg1Pos2 || B == LUFpNormalizeBound::HiNeg4Pos4);
        }
    }

    static_assert(std::is_trivially_default_constructible_v<FpPack8x1>);
    static_assert(std::is_trivially_copyable_v<FpPack8x1>);
    static_assert(std::is_trivially_default_constructible_v<LFpPack8x1>);
    static_assert(std::is_trivially_copyable_v<LFpPack8x1>);
    static_assert(std::is_trivially_default_constructible_v<LUFpPack8x1>);
    static_assert(std::is_trivially_copyable_v<LUFpPack8x1>);
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


    /// @brief Generic multi-precision conversion from base 2^64 to base 2^48.
    /// @param rlen Number of output radix-48 limbs.
    /// @param alen Number of input 64-bit limbs.
    /// @post High unused output limbs are zero-filled.
    void conv_64to48_mpi(u64* r, const u64* a, int rlen, int alen);


    /// @brief Convert one scalar from 8x48-bit limbs to 6x64-bit limbs.
    /// @pre `a` points to 8 radix-48 limbs.
    /// @post `r` reconstructs the same integer in 64-bit limb form.
    void conv_48to64_fp(u64* r, const u64* a);


    /// @brief Generic multi-precision conversion from base 2^48 to base 2^64.
    /// @post High unused output limbs are zero-filled.
    void conv_48to64_mpi(u64* r, const u64* a, int rlen, int alen);


    /// @brief Carry-propagate scalar radix-48 limbs in place.
    /// @pre `a` has at least `NWORDS` limbs.
    /// @post `a[0..NWORDS-2] < 2^48`; final carry remains in `a[NWORDS-1]`.
    void carryp_mpi48(u64* a);
    ///@}

    /// @brief Debug-only telemetry for observed lazy-reduced limb bounds.
    struct LFpBoundsStats8x1
    {
        u64 total_lanes;
        u64 ge_p_lanes;
        u64 max_limb;
    };

    /// @brief Reset lazy-bound telemetry counters.
    void fp8x1_lfp_bounds_reset();

    /// @brief Read lazy-bound telemetry counters.
    LFpBoundsStats8x1 fp8x1_lfp_bounds_get();

    /// @brief K-wise widened multiply over `FpPack8x1` arrays.
    /// @tparam K Number of independent packs processed per block (`1`, `2`, or `4`).
    /// @param n_blocks Number of K-sized blocks to process.
    template <int K>
    BLST_AVX_NOINLINE void fp8x1_mul_k(
        UFpPack8x1* BLST_restrict out,
        const FpPack8x1* BLST_restrict a,
        const FpPack8x1* BLST_restrict b,
        std::size_t n_blocks);

    /// @brief K-wise Montgomery reduction over widened `UFpPack8x1` arrays.
    /// @tparam K Number of independent packs processed per block (`1`, `2`, or `4`).
    /// @param n_blocks Number of K-sized blocks to process.
    template <int K>
    BLST_AVX_NOINLINE void fp8x1_redc_k(
        FpPack8x1* BLST_restrict out,
        const UFpPack8x1* BLST_restrict in,
        std::size_t n_blocks);

    /// @brief K-wise composed multiply+reduce over `FpPack8x1` arrays.
    /// @tparam K Number of independent packs processed per block (`1`, `2`, or `4`).
    /// @param n_blocks Number of K-sized blocks to process.
    template <int K>
    BLST_AVX_NOINLINE void fp8x1_mul_redc_k(
        FpPack8x1* BLST_restrict out,
        const FpPack8x1* BLST_restrict a,
        const FpPack8x1* BLST_restrict b,
        std::size_t n_blocks);

    /// @brief Runtime 4->2->1 dispatch for widened multiply.
    void fp8x1_mul_dispatch(
        UFpPack8x1* BLST_restrict out,
        const FpPack8x1* BLST_restrict a,
        const FpPack8x1* BLST_restrict b,
        std::size_t n_packs);

    /// @brief Runtime 4->2->1 dispatch for Montgomery reduction.
    void fp8x1_redc_dispatch(
        FpPack8x1* BLST_restrict out,
        const UFpPack8x1* BLST_restrict in,
        std::size_t n_packs);

    /// @brief Runtime 4->2->1 dispatch for composed multiply+reduce.
    void fp8x1_mul_redc_dispatch(
        FpPack8x1* BLST_restrict out,
        const FpPack8x1* BLST_restrict a,
        const FpPack8x1* BLST_restrict b,
        std::size_t n_packs);

    BLST_FORCEINLINE void LUFpPack8x1::add(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b)
    {
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            r[i] = VADD(a[i], b[i]);
        }
    }

    BLST_FORCEINLINE void LUFpPack8x1::add(LUFpPack8x1& r, const UFpPack8x1& a, const LUFpPack8x1& b)
    {
        LUFpPack8x1::add(r, reinterpret_cast<const LUFpPack8x1&>(a), b);
    }

    BLST_FORCEINLINE void LUFpPack8x1::add(LUFpPack8x1& r, const LUFpPack8x1& a, const UFpPack8x1& b)
    {
        LUFpPack8x1::add(r, a, reinterpret_cast<const LUFpPack8x1&>(b));
    }

    BLST_FORCEINLINE void LUFpPack8x1::add(LUFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        LUFpPack8x1::add(r, reinterpret_cast<const LUFpPack8x1&>(a), reinterpret_cast<const LUFpPack8x1&>(b));
    }

    BLST_FORCEINLINE void LUFpPack8x1::sub(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b)
    {
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            r[i] = VSUB(a[i], b[i]);
        }
    }

    BLST_FORCEINLINE void LUFpPack8x1::sub(LUFpPack8x1& r, const UFpPack8x1& a, const LUFpPack8x1& b)
    {
        LUFpPack8x1::sub(r, reinterpret_cast<const LUFpPack8x1&>(a), b);
    }

    BLST_FORCEINLINE void LUFpPack8x1::sub(LUFpPack8x1& r, const LUFpPack8x1& a, const UFpPack8x1& b)
    {
        LUFpPack8x1::sub(r, a, reinterpret_cast<const LUFpPack8x1&>(b));
    }

    BLST_FORCEINLINE void LUFpPack8x1::sub(LUFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        LUFpPack8x1::sub(r, reinterpret_cast<const LUFpPack8x1&>(a), reinterpret_cast<const LUFpPack8x1&>(b));
    }

    inline void UFpPack8x1::add(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
        const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
        const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
        const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
        const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
        const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
        const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
        const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
        cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask, t0;

        r0  = VADD(a0 , b0 ); r1  = VADD(a1 , b1 );
        r2  = VADD(a2 , b2 ); r3  = VADD(a3 , b3 );
        r4  = VADD(a4 , b4 ); r5  = VADD(a5 , b5 );
        r6  = VADD(a6 , b6 ); r7  = VADD(a7 , b7 );
        r8  = VADD(a8 , b8 ); r9  = VADD(a9 , b9 );
        r10 = VADD(a10, b10); r11 = VADD(a11, b11);
        r12 = VADD(a12, b12); r13 = VADD(a13, b13);
        r14 = VADD(a14, b14); r15 = VADD(a15, b15);

        r8  = VSUB(r8 , p0); r9  = VSUB(r9 , p1);
        r10 = VSUB(r10, p2); r11 = VSUB(r11, p3);
        r12 = VSUB(r12, p4); r13 = VSUB(r13, p5);
        r14 = VSUB(r14, p6); r15 = VSUB(r15, p7);

        r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
        t0 = VADD(r9 , VSRA(r8, BRADIX));
        t0 = VADD(r10, VSRA(t0, BRADIX));
        t0 = VADD(r11, VSRA(t0, BRADIX));
        t0 = VADD(r12, VSRA(t0, BRADIX));
        t0 = VADD(r13, VSRA(t0, BRADIX));
        t0 = VADD(r14, VSRA(t0, BRADIX));
        t0 = VADD(r15, VSRA(t0, BRADIX));

        smask = VSRA(t0, 63);
        r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
        r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
        r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
        r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

        r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
        r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
        r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
        r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
        r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
        r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
        r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                          r15 = VAND(r15, bmask);

        r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
        r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
        r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
        r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
    }

    inline void UFpPack8x1::sub(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
        const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
        const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
        const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
        const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
        const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
        const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
        const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
        cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask, t0;

        r0  = VSUB(a0 , b0 ); r1  = VSUB(a1 , b1 );
        r2  = VSUB(a2 , b2 ); r3  = VSUB(a3 , b3 );
        r4  = VSUB(a4 , b4 ); r5  = VSUB(a5 , b5 );
        r6  = VSUB(a6 , b6 ); r7  = VSUB(a7 , b7 );
        r8  = VSUB(a8 , b8 ); r9  = VSUB(a9 , b9 );
        r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
        r12 = VSUB(a12, b12); r13 = VSUB(a13, b13);
        r14 = VSUB(a14, b14); r15 = VSUB(a15, b15);

        r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
        t0 = VADD(r9 , VSRA(r8, BRADIX));
        t0 = VADD(r10, VSRA(t0, BRADIX));
        t0 = VADD(r11, VSRA(t0, BRADIX));
        t0 = VADD(r12, VSRA(t0, BRADIX));
        t0 = VADD(r13, VSRA(t0, BRADIX));
        t0 = VADD(r14, VSRA(t0, BRADIX));
        t0 = VADD(r15, VSRA(t0, BRADIX));

        smask = VSRA(t0, 63);
        r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
        r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
        r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
        r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

        r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
        r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
        r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
        r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
        r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
        r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
        r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                          r15 = VAND(r15, bmask);

        r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
        r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
        r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
        r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
    }


    inline void FpPack8x1::add(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask, t0;

        // r = a + b
        r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); r3 = VADD(a3, b3);
        r4 = VADD(a4, b4); r5 = VADD(a5, b5); r6 = VADD(a6, b6); r7 = VADD(a7, b7);

        // r = r - p
        r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);
        r4 = VSUB(r4, p4); r5 = VSUB(r5, p5); r6 = VSUB(r6, p6); r7 = VSUB(r7, p7);

        // get sign mask 
        t0 = VADD(r1, VSRA(r0, BRADIX));
        t0 = VADD(r2, VSRA(t0, BRADIX));
        t0 = VADD(r3, VSRA(t0, BRADIX));
        t0 = VADD(r4, VSRA(t0, BRADIX));
        t0 = VADD(r5, VSRA(t0, BRADIX));
        t0 = VADD(r6, VSRA(t0, BRADIX));
        t0 = VADD(r7, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
        r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
        r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

        // carry propagation 
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r7 = VAND(r7, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
        r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    }

    inline void FpPack8x1::sub(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask, t0;

        // r = a - b
        r0 = VSUB(a0, b0); r1 = VSUB(a1, b1); r2 = VSUB(a2, b2); r3 = VSUB(a3, b3);
        r4 = VSUB(a4, b4); r5 = VSUB(a5, b5); r6 = VSUB(a6, b6); r7 = VSUB(a7, b7);

        // get sign mask 
        t0 = VADD(r1, VSRA(r0, BRADIX));
        t0 = VADD(r2, VSRA(t0, BRADIX));
        t0 = VADD(r3, VSRA(t0, BRADIX));
        t0 = VADD(r4, VSRA(t0, BRADIX));
        t0 = VADD(r5, VSRA(t0, BRADIX));
        t0 = VADD(r6, VSRA(t0, BRADIX));
        t0 = VADD(r7, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
        r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
        r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

        // carry propagation 
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r7 = VAND(r7, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
        r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    }



    inline void fp8x1_mul_schoolbook_paircarry(
        UFpPack8x1& r,
        const FpPack8x1& a,
        const FpPack8x1& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];

        cc512i carry = VZERO();
        cc512i lo0 = VZERO(), hi0 = VZERO();
        cc512i lo1 = VZERO(), hi1 = VZERO();

#define BLST_SB_PC_PAIR_BEGIN() \
  lo0 = carry;                  \
  hi0 = VZERO()
#define BLST_SB_PC_TERM0(aa, bb)    \
  lo0 = VMACLO(lo0, (aa), (bb));    \
  hi0 = VMACHI(hi0, (aa), (bb))
#define BLST_SB_PC_PAIR_MID() \
  lo1 = VSHL(hi0, BALIGN);    \
  hi1 = VZERO()
#define BLST_SB_PC_TERM1(aa, bb)    \
  lo1 = VMACLO(lo1, (aa), (bb));    \
  hi1 = VMACHI(hi1, (aa), (bb))
#define BLST_SB_PC_PAIR_END(idx0, idx1) \
  r[(idx0)] = lo0;                      \
  r[(idx1)] = lo1;                      \
  carry = VSHL(hi1, BALIGN)

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a0, b0);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a0, b1); BLST_SB_PC_TERM1(a1, b0);
        BLST_SB_PC_PAIR_END(0, 1);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a0, b2); BLST_SB_PC_TERM0(a1, b1); BLST_SB_PC_TERM0(a2, b0);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a0, b3); BLST_SB_PC_TERM1(a1, b2); BLST_SB_PC_TERM1(a2, b1); BLST_SB_PC_TERM1(a3, b0);
        BLST_SB_PC_PAIR_END(2, 3);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a0, b4); BLST_SB_PC_TERM0(a1, b3); BLST_SB_PC_TERM0(a2, b2); BLST_SB_PC_TERM0(a3, b1);
        BLST_SB_PC_TERM0(a4, b0);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a0, b5); BLST_SB_PC_TERM1(a1, b4); BLST_SB_PC_TERM1(a2, b3); BLST_SB_PC_TERM1(a3, b2);
        BLST_SB_PC_TERM1(a4, b1); BLST_SB_PC_TERM1(a5, b0);
        BLST_SB_PC_PAIR_END(4, 5);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a0, b6); BLST_SB_PC_TERM0(a1, b5); BLST_SB_PC_TERM0(a2, b4); BLST_SB_PC_TERM0(a3, b3);
        BLST_SB_PC_TERM0(a4, b2); BLST_SB_PC_TERM0(a5, b1); BLST_SB_PC_TERM0(a6, b0);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a0, b7); BLST_SB_PC_TERM1(a1, b6); BLST_SB_PC_TERM1(a2, b5); BLST_SB_PC_TERM1(a3, b4);
        BLST_SB_PC_TERM1(a4, b3); BLST_SB_PC_TERM1(a5, b2); BLST_SB_PC_TERM1(a6, b1); BLST_SB_PC_TERM1(a7, b0);
        BLST_SB_PC_PAIR_END(6, 7);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a1, b7); BLST_SB_PC_TERM0(a2, b6); BLST_SB_PC_TERM0(a3, b5); BLST_SB_PC_TERM0(a4, b4);
        BLST_SB_PC_TERM0(a5, b3); BLST_SB_PC_TERM0(a6, b2); BLST_SB_PC_TERM0(a7, b1);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a2, b7); BLST_SB_PC_TERM1(a3, b6); BLST_SB_PC_TERM1(a4, b5); BLST_SB_PC_TERM1(a5, b4);
        BLST_SB_PC_TERM1(a6, b3); BLST_SB_PC_TERM1(a7, b2);
        BLST_SB_PC_PAIR_END(8, 9);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a3, b7); BLST_SB_PC_TERM0(a4, b6); BLST_SB_PC_TERM0(a5, b5); BLST_SB_PC_TERM0(a6, b4);
        BLST_SB_PC_TERM0(a7, b3);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a4, b7); BLST_SB_PC_TERM1(a5, b6); BLST_SB_PC_TERM1(a6, b5); BLST_SB_PC_TERM1(a7, b4);
        BLST_SB_PC_PAIR_END(10, 11);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a5, b7); BLST_SB_PC_TERM0(a6, b6); BLST_SB_PC_TERM0(a7, b5);
        BLST_SB_PC_PAIR_MID();
        BLST_SB_PC_TERM1(a6, b7); BLST_SB_PC_TERM1(a7, b6);
        BLST_SB_PC_PAIR_END(12, 13);

        BLST_SB_PC_PAIR_BEGIN();
        BLST_SB_PC_TERM0(a7, b7);
        r[14] = lo0;
        r[15] = VSHL(hi0, BALIGN);

#undef BLST_SB_PC_PAIR_BEGIN
#undef BLST_SB_PC_TERM0
#undef BLST_SB_PC_PAIR_MID
#undef BLST_SB_PC_TERM1
#undef BLST_SB_PC_PAIR_END
    }

    BLST_FORCEINLINE void FpPack8x1::mul(UFpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b)
    {
        BLST_ASSERT_NO_OVERLAP_AVX(r, a);
        BLST_ASSERT_NO_OVERLAP_AVX(r, b);
        fp8x1_mul_schoolbook_paircarry(r, a, b);
    }


}

namespace blst
{
    using avx::BRADIX;
    using avx::BALIGN;
    using avx::BMASK;
    using avx::MONT_W_R48;
    using avx::NWORDS;
    using avx::VWORDS;
    using avx::conv_64to48_fp;
    using avx::conv_48to64_fp;
    using avx::carryp_mpi48;
    using avx::LFpBoundsStats8x1;
    using avx::fp8x1_lfp_bounds_reset;
    using avx::fp8x1_lfp_bounds_get;
    using avx::fp8x1_mul_dispatch;
    using avx::fp8x1_redc_dispatch;
    using avx::fp8x1_mul_redc_dispatch;

    using FpPack8 = avx::FpPack8x1;
    using LFpPack8 = avx::LFpPack8x1;
    using UFpPack8 = avx::UFpPack8x1;
    using LUFpPack8 = avx::LUFpPack8x1;

    using FpPack8x1 = avx::FpPack8x1;
    using LFpPack8x1 = avx::LFpPack8x1;
    using UFpPack8x1 = avx::UFpPack8x1;
    using LUFpPack8x1 = avx::LUFpPack8x1;
    using FpPack4x2 = avx::FpPack4x2;
    using UFpPack4x2 = avx::UFpPack4x2;
    using LUFpRenormClass = avx::LUFpRenormClass;
    using LUFpNormalizeBound = avx::LUFpNormalizeBound;
    extern const std::array<u64, avx::NWORDS> P48;

    template<int K>
    BLST_AVX_NOINLINE void fp8x1_mul_k(cc512i* BLST_restrict r, const cc512i* BLST_restrict a,
        const cc512i* BLST_restrict b, std::size_t n_blocks)
    {
        avx::fp8x1_mul_k<K>(r, a, b, n_blocks);
    }

    template<int K>
    BLST_AVX_NOINLINE void fp8x1_redc_k(cc512i* BLST_restrict r, const cc512i* BLST_restrict a,
        std::size_t n_blocks)
    {
        avx::fp8x1_redc_k<K>(r, a, n_blocks);
    }

    template<int K>
    BLST_AVX_NOINLINE void fp8x1_mul_redc_k(cc512i* BLST_restrict r, const cc512i* BLST_restrict a,
        const cc512i* BLST_restrict b, std::size_t n_blocks)
    {
        avx::fp8x1_mul_redc_k<K>(r, a, b, n_blocks);
    }
}

#endif

#endif // BLST_FIELDS_FP_AVX_H
