#ifndef BLST_FIELDS_FP_AVX_H
#define BLST_FIELDS_FP_AVX_H

#include "blst/common/defines.h"
#include "blst/backend/avx512/avx_ops.h"
#include <cstddef>
#include <type_traits>

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

namespace blst
{
    struct Fp;

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
    };


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
        } else {
            static_assert(B == LUFpNormalizeBound::HiNeg1Pos2);
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
}

#endif

#endif // BLST_FIELDS_FP_AVX_H
