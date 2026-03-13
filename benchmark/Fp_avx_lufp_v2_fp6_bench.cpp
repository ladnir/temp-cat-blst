#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "blst/fields/Fp.h"
#include "benchmark/internal/Fp_avx_lufp_v2_impl.h"

using namespace blst;
using namespace blst::bench_lufp_v2;

namespace
{
    bool fp2_equal(const Fp2& a, const Fp2& b)
    {
        return Fp::equal(a[0], b[0]) && Fp::equal(a[1], b[1]);
    }

    bool fp6_equal(const Fp6& a, const Fp6& b)
    {
        return fp2_equal(a[0], b[0]) && fp2_equal(a[1], b[1]) && fp2_equal(a[2], b[2]);
    }

    bool fp6_wide_equal(const UFp6Pack8x1& a, const UFp6Pack8x1& b)
    {
        for (std::size_t c = 0; c < 3; ++c) {
            for (std::size_t d = 0; d < 2; ++d) {
                if (!lufp_v2_equal(a.c[c][d], b.c[c][d])) {
                    return false;
                }
            }
        }
        return true;
    }

    void dump_fp(const char* label, const Fp& v)
    {
        std::fprintf(
            stderr,
            "%s = [%016llx %016llx %016llx %016llx %016llx %016llx]\n",
            label,
            static_cast<unsigned long long>(v[5]),
            static_cast<unsigned long long>(v[4]),
            static_cast<unsigned long long>(v[3]),
            static_cast<unsigned long long>(v[2]),
            static_cast<unsigned long long>(v[1]),
            static_cast<unsigned long long>(v[0]));
    }

    void dump_fp2(const char* label, const Fp2& v)
    {
        std::fprintf(stderr, "%s.re\n", label);
        dump_fp("  ", v[0]);
        std::fprintf(stderr, "%s.im\n", label);
        dump_fp("  ", v[1]);
    }

    void dump_fp6(const char* label, const Fp6& v)
    {
        std::fprintf(stderr, "%s.c0\n", label);
        dump_fp2("  ", v[0]);
        std::fprintf(stderr, "%s.c1\n", label);
        dump_fp2("  ", v[1]);
        std::fprintf(stderr, "%s.c2\n", label);
        dump_fp2("  ", v[2]);
    }

    using PackedI64Words = std::array<std::array<i64, 8>, 2 * NWORDS>;

    template <class Pack>
    PackedI64Words pack_to_i64_words(const Pack& in)
    {
        PackedI64Words out{};
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            std::array<u64, 8> lane_words{};
            _cc512_storeu_si512(lane_words.data(), in[limb]);
            for (std::size_t lane = 0; lane < lane_words.size(); ++lane) {
                out[limb][lane] = static_cast<i64>(lane_words[lane]);
            }
        }
        return out;
    }

    template <class Pack>
    void pack_from_i64_words(Pack& out, const PackedI64Words& in)
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            std::array<u64, 8> lane_words{};
            for (std::size_t lane = 0; lane < lane_words.size(); ++lane) {
                lane_words[lane] = static_cast<u64>(in[limb][lane]);
            }
            out[limb] = _cc512_loadu_si512(lane_words.data());
        }
    }

    bool packed_words_equal(const PackedI64Words& a, const PackedI64Words& b)
    {
        for (std::size_t limb = 0; limb < a.size(); ++limb) {
            for (std::size_t lane = 0; lane < a[limb].size(); ++lane) {
                if (a[limb][lane] != b[limb][lane]) {
                    return false;
                }
            }
        }
        return true;
    }

    void dump_packed_words(const char* label, const PackedI64Words& words)
    {
        std::fprintf(stderr, "%s\n", label);
        for (std::size_t limb = 0; limb < words.size(); ++limb) {
            std::fprintf(stderr, "  limb[%zu] =", limb);
            for (std::size_t lane = 0; lane < words[limb].size(); ++lane) {
                std::fprintf(stderr, " %lld", static_cast<long long>(words[limb][lane]));
            }
            std::fprintf(stderr, "\n");
        }
    }

    u8 high_negative_mask_scalar(const PackedI64Words& words)
    {
        u8 mask = 0;
        for (u8 lane = 0; lane < 8; ++lane) {
            mask |= static_cast<u8>(((words[(2 * NWORDS) - 1][lane] < 0) ? 1u : 0u) << lane);
        }
        return mask;
    }

    v2_lane_mask_t u8_to_lane_mask(const u8 mask)
    {
#if BLST_AVX_NATIVE
        return static_cast<v2_lane_mask_t>(mask);
#else
        return mask;
#endif
    }

    u8 lane_mask_to_u8(const v2_lane_mask_t mask)
    {
#if BLST_AVX_NATIVE
        return static_cast<u8>(mask);
#else
        return mask;
#endif
    }

    i64 scalar_floor_div_radix(const i64 x)
    {
        if (x >= 0) {
            return x >> BRADIX;
        }
        return -static_cast<i64>((static_cast<u64>(-x) + BMASK) >> BRADIX);
    }

    void scalar_signed_carry_range(PackedI64Words& words, const std::size_t begin, const std::size_t end)
    {
        for (std::size_t limb = begin; limb < end; ++limb) {
            for (std::size_t lane = 0; lane < 8; ++lane) {
                const i64 carry = scalar_floor_div_radix(words[limb][lane]);
                words[limb][lane] -= (carry << BRADIX);
                words[limb + 1][lane] += carry;
            }
        }
    }

    void dump_step_candidate_half_debug(const char* label, const UFpPack8x1& raw)
    {
        UFpPack8x1 safe_state = raw;
        UFpPack8x1 candidate_state = raw;

        v2_signed_carry_range<0, (2 * NWORDS) - 1>(safe_state.data());
        const u8 neg_mask = lane_mask_to_u8(v2_signbit_mask(safe_state[(2 * NWORDS) - 1]));
        v2_high_add_mod_masked_u8(safe_state.data(), neg_mask, kV2P48X4.data());
        v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(safe_state.data());
        const u8 neg_fix_mask = lane_mask_to_u8(v2_signbit_mask(safe_state[(2 * NWORDS) - 1]));
        v2_high_add_mod_masked_u8(safe_state.data(), neg_fix_mask, kV2P48X1.data());
        v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(safe_state.data());
        const u8 ge_4p_mask = v2_high_ge_mod_mask_u8(safe_state.data(), kV2P48X4.data());
        const u8 ge_2p_mask = v2_high_ge_mod_mask_u8(safe_state.data(), kV2P48X2.data());
        const u8 ge_p_mask = v2_high_ge_p_mask_u8(safe_state.data());

        v2_high_sub_const(candidate_state.data(), kV2Bias3PPlus3Vec);
        const u8 sign0 = lane_mask_to_u8(v2_carry_range_with_sign<0, (2 * NWORDS) - 1>(candidate_state.data()));
        const u8 sign4 = lane_mask_to_u8(v2_sign_step<4>(candidate_state.data(), static_cast<v2_lane_mask_t>(sign0)));
        const u8 sign42 = lane_mask_to_u8(v2_sign_step<2>(candidate_state.data(), static_cast<v2_lane_mask_t>(sign4)));
        const u8 sign421 = lane_mask_to_u8(v2_sign_step<1>(candidate_state.data(), static_cast<v2_lane_mask_t>(sign42)));
        const u8 sign4211 = lane_mask_to_u8(v2_sign_step<1>(candidate_state.data(), static_cast<v2_lane_mask_t>(sign421)));

        std::fprintf(
            stderr,
            "%s masks lane0: safe neg=%u neg_fix=%u ge4=%u ge2=%u gep=%u | cand s0=%u s4=%u s42=%u s421=%u s4211=%u\n",
            label,
            static_cast<unsigned>(neg_mask & 1u),
            static_cast<unsigned>(neg_fix_mask & 1u),
            static_cast<unsigned>(ge_4p_mask & 1u),
            static_cast<unsigned>(ge_2p_mask & 1u),
            static_cast<unsigned>(ge_p_mask & 1u),
            static_cast<unsigned>(sign0 & 1u),
            static_cast<unsigned>(sign4 & 1u),
            static_cast<unsigned>(sign42 & 1u),
            static_cast<unsigned>(sign421 & 1u),
            static_cast<unsigned>(sign4211 & 1u));
    }

    void assert_fp2_pair_equal_or_dump(const char* label, const UFp2Pack8x1x1& got, const UFp2Pack8x1x1& want)
    {
        const PackedI64Words got_re = pack_to_i64_words(got[0]);
        const PackedI64Words want_re = pack_to_i64_words(want[0]);
        const PackedI64Words got_im = pack_to_i64_words(got[1]);
        const PackedI64Words want_im = pack_to_i64_words(want[1]);
        if (!packed_words_equal(got_re, want_re) || !packed_words_equal(got_im, want_im)) {
            std::fprintf(stderr, "step divergence at %s\n", label);
            dump_packed_words("got.re", got_re);
            dump_packed_words("want.re", want_re);
            dump_packed_words("got.im", got_im);
            dump_packed_words("want.im", want_im);
            std::abort();
        }
    }

    bool fp2_pair_words_equal(const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        return packed_words_equal(pack_to_i64_words(a[0]), pack_to_i64_words(b[0])) &&
            packed_words_equal(pack_to_i64_words(a[1]), pack_to_i64_words(b[1]));
    }

    void report_half_repair_match(const char* half_label, const char* repair_label, const UFpPack8x1& got, const UFpPack8x1& want)
    {
        if (packed_words_equal(pack_to_i64_words(got), pack_to_i64_words(want))) {
            std::fprintf(stderr, "  %s matches %s\n", half_label, repair_label);
        }
    }

    void debug_search_half_sign_sequences(const char* half_label, const UFpPack8x1& raw, const UFpPack8x1& want)
    {
        static constexpr u64 kScales[3] = { 1, 2, 4 };
        static constexpr const char* kScaleNames[3] = { "1", "2", "4" };

        auto print_sequence = [&](const std::array<int, 6>& digits, const int len) {
            std::fprintf(stderr, "  %s sign-seq:", half_label);
            if (len == 0) {
                std::fprintf(stderr, " <empty>");
            }
            for (int i = 0; i < len; ++i) {
                std::fprintf(stderr, " %s", kScaleNames[digits[static_cast<std::size_t>(i)]]);
            }
            std::fprintf(stderr, "\n");
            };

        for (int len = 0; len <= 6; ++len) {
            const int variants = len == 0 ? 1 : static_cast<int>(std::pow(3, len));
            int printed = 0;
            for (int variant = 0; variant < variants; ++variant) {
                std::array<int, 6> digits{};
                int code = variant;
                for (int i = 0; i < len; ++i) {
                    digits[static_cast<std::size_t>(i)] = code % 3;
                    code /= 3;
                }

                UFpPack8x1 got = raw;
                v2_high_sub_const(got.data(), kV2Bias3PPlus3Vec);
                v2_signed_carry_range<0, (2 * NWORDS) - 1>(got.data());
                v2_lane_mask_t sign_mask = v2_signbit_mask(got[(2 * NWORDS) - 1]);
                for (int i = 0; i < len; ++i) {
                    const int digit = digits[static_cast<std::size_t>(i)];
                    if (kScales[digit] == 1) {
                        sign_mask = v2_sign_step<1>(got.data(), sign_mask);
                    }
                    else if (kScales[digit] == 2) {
                        sign_mask = v2_sign_step<2>(got.data(), sign_mask);
                    }
                    else {
                        sign_mask = v2_sign_step<4>(got.data(), sign_mask);
                    }
                }

                if (packed_words_equal(pack_to_i64_words(got), pack_to_i64_words(want))) {
                    print_sequence(digits, len);
                    ++printed;
                    if (printed == 4) {
                        break;
                    }
                }
            }
            if (printed != 0) {
                return;
            }
        }
        std::fprintf(stderr, "  %s sign-seq: no match up to len 6\n", half_label);
    }

    void lufp2_v2_renorm_exact_add_pair(UFp2Pack8x1x1& in_out)
    {
        lufp_v2_renorm_exact_add(in_out[0]);
        lufp_v2_renorm_exact_add(in_out[1]);
    }

    u32 popcount_u32(u32 value)
    {
        u32 count = 0;
        while (value != 0) {
            count += (value & 1u);
            value >>= 1;
        }
        return count;
    }

    void maybe_renorm_exact_add_halves(UFp2Pack8x1x1& in_out, const u32 re_mask, const u32 im_mask, const u32 bit)
    {
        if ((re_mask & (1u << bit)) != 0u) {
            lufp_v2_renorm_exact_add(in_out[0]);
        }
        if ((im_mask & (1u << bit)) != 0u) {
            lufp_v2_renorm_exact_add(in_out[1]);
        }
    }

    bool step_candidate_half_mask_matches_reference(const Fp6Pack8x1& a, const Fp6Pack8x1& b, const u32 re_mask, const u32 im_mask)
    {
        Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
        Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
        UFp2Pack8x1x1 t0{}, t1{}, t2{};
        UFp2Pack8x1x1 w0{}, w1{}, w2{};
        UFp2Pack8x1x1 t3{};
        UFp6Pack8x1 candidate{};
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 ref_mul{};
        std::array<Fp6, 8> got_candidate{};
        std::array<Fp6, 8> ref{};

        Fp2Pack8x1x1::add(aa0, a.c[1], a.c[2]);
        Fp2Pack8x1x1::add(bb0, b.c[1], b.c[2]);
        Fp2Pack8x1x1::add(aa1, a.c[0], a.c[1]);
        Fp2Pack8x1x1::add(bb1, b.c[0], b.c[1]);
        Fp2Pack8x1x1::add(aa2, a.c[0], a.c[2]);
        Fp2Pack8x1x1::add(bb2, b.c[0], b.c[2]);

        UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
        UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
        UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);
        UFp2Pack8x1x1::mul(w0, aa0, bb0);
        UFp2Pack8x1x1::mul(w1, aa1, bb1);
        UFp2Pack8x1x1::mul(w2, aa2, bb2);

        candidate.c[0] = w0;
        lufp2_v2_sub(candidate.c[0], candidate.c[0], t1);
        maybe_renorm_exact_add_halves(candidate.c[0], re_mask, im_mask, 0);
        lufp2_v2_sub(candidate.c[0], candidate.c[0], t2);
        maybe_renorm_exact_add_halves(candidate.c[0], re_mask, im_mask, 1);
        lufp2_v2_mul_by_u_plus_1(t3, candidate.c[0]);
        maybe_renorm_exact_add_halves(t3, re_mask, im_mask, 2);
        lufp2_v2_add(candidate.c[0], t3, t0);
        maybe_renorm_exact_add_halves(candidate.c[0], re_mask, im_mask, 3);

        candidate.c[1] = w1;
        lufp2_v2_sub(candidate.c[1], candidate.c[1], t0);
        maybe_renorm_exact_add_halves(candidate.c[1], re_mask, im_mask, 4);
        lufp2_v2_sub(candidate.c[1], candidate.c[1], t1);
        maybe_renorm_exact_add_halves(candidate.c[1], re_mask, im_mask, 5);
        lufp2_v2_mul_by_u_plus_1(t3, t2);
        maybe_renorm_exact_add_halves(t3, re_mask, im_mask, 6);
        lufp2_v2_add(candidate.c[1], candidate.c[1], t3);
        maybe_renorm_exact_add_halves(candidate.c[1], re_mask, im_mask, 7);

        candidate.c[2] = w2;
        lufp2_v2_sub(candidate.c[2], candidate.c[2], t0);
        maybe_renorm_exact_add_halves(candidate.c[2], re_mask, im_mask, 8);
        lufp2_v2_sub(candidate.c[2], candidate.c[2], t2);
        maybe_renorm_exact_add_halves(candidate.c[2], re_mask, im_mask, 9);
        lufp2_v2_add(candidate.c[2], candidate.c[2], t1);
        maybe_renorm_exact_add_halves(candidate.c[2], re_mask, im_mask, 10);

        UFp6Pack8x1::redc(candidate_redc, candidate);
        Fp6Pack8x1::mul(ref_mul, a, b);
        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);
        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                return false;
            }
        }
        return true;
    }

    void debug_search_step_candidate_masks(const Fp6Pack8x1& a, const Fp6Pack8x1& b)
    {
        static constexpr const char* kCheckpointNames[11] = {
            "c0_sub_t1",
            "c0_sub_t2",
            "c0_uplus1",
            "c0_add_t0",
            "c1_sub_t0",
            "c1_sub_t1",
            "c1_uplus1_t2",
            "c1_add_t3",
            "c2_sub_t0",
            "c2_sub_t2",
            "c2_add_t1",
        };

        u32 best_count = 32;
        std::array<u32, 16> best_masks{};
        std::size_t best_masks_size = 0;

        for (u32 mask = 0; mask < (1u << 11); ++mask) {
            if (!step_candidate_half_mask_matches_reference(a, b, mask, mask)) {
                continue;
            }

            const u32 count = popcount_u32(mask);
            if (count < best_count) {
                best_count = count;
                best_masks_size = 0;
            }
            if (count == best_count && best_masks_size < best_masks.size()) {
                best_masks[best_masks_size++] = mask;
            }
        }

        if (best_masks_size == 0) {
            std::fprintf(stderr, "no passing exact-add subset found\n");
            return;
        }

        std::fprintf(stderr, "best passing exact-add subset count: %u\n", static_cast<unsigned>(best_count));
        for (std::size_t i = 0; i < best_masks_size; ++i) {
            const u32 mask = best_masks[i];
            std::fprintf(stderr, "  mask 0x%03x:", static_cast<unsigned>(mask));
            for (u32 bit = 0; bit < 11; ++bit) {
                if ((mask & (1u << bit)) != 0u) {
                    std::fprintf(stderr, " %s", kCheckpointNames[bit]);
                }
            }
            std::fprintf(stderr, "\n");
        }
    }

    void debug_search_step_candidate_half_masks(const Fp6Pack8x1& a, const Fp6Pack8x1& b)
    {
        static constexpr const char* kCheckpointNames[11] = {
            "c0_sub_t1",
            "c0_sub_t2",
            "c0_uplus1",
            "c0_add_t0",
            "c1_sub_t0",
            "c1_sub_t1",
            "c1_uplus1_t2",
            "c1_add_t3",
            "c2_sub_t0",
            "c2_sub_t2",
            "c2_add_t1",
        };
        static constexpr u32 kBaseMasks[4] = { 0x362u, 0x3C2u, 0x562u, 0x5C2u };

        u32 best_count = 32;
        struct HalfMaskResult
        {
            u32 base_mask;
            u32 re_mask;
            u32 im_mask;
        };
        std::array<HalfMaskResult, 16> best_results{};
        std::size_t best_results_size = 0;

        for (const u32 base_mask : kBaseMasks) {
            for (u32 re_mask = base_mask;; re_mask = (re_mask - 1u) & base_mask) {
                for (u32 im_mask = base_mask;; im_mask = (im_mask - 1u) & base_mask) {
                    if (step_candidate_half_mask_matches_reference(a, b, re_mask, im_mask)) {
                        const u32 count = popcount_u32(re_mask) + popcount_u32(im_mask);
                        if (count < best_count) {
                            best_count = count;
                            best_results_size = 0;
                        }
                        if (count == best_count && best_results_size < best_results.size()) {
                            best_results[best_results_size++] = { base_mask, re_mask, im_mask };
                        }
                    }
                    if (im_mask == 0u) {
                        break;
                    }
                }
                if (re_mask == 0u) {
                    break;
                }
            }
        }

        if (best_results_size == 0) {
            std::fprintf(stderr, "no passing half exact-add subset found under sparse-5 masks\n");
            return;
        }

        std::fprintf(stderr, "best passing half exact-add count under sparse-5 masks: %u\n", static_cast<unsigned>(best_count));
        for (std::size_t i = 0; i < best_results_size; ++i) {
            const HalfMaskResult& result = best_results[i];
            std::fprintf(stderr, "  base 0x%03x | re:", static_cast<unsigned>(result.base_mask));
            for (u32 bit = 0; bit < 11; ++bit) {
                if ((result.re_mask & (1u << bit)) != 0u) {
                    std::fprintf(stderr, " %s", kCheckpointNames[bit]);
                }
            }
            std::fprintf(stderr, " | im:");
            for (u32 bit = 0; bit < 11; ++bit) {
                if ((result.im_mask & (1u << bit)) != 0u) {
                    std::fprintf(stderr, " %s", kCheckpointNames[bit]);
                }
            }
            std::fprintf(stderr, "\n");
        }
    }

    template <u64 Scale>
    void validate_v2_sign_step_case(std::mt19937_64& gen)
    {
        constexpr i64 kRange = static_cast<i64>(1ULL << 47);
        std::uniform_int_distribution<i64> dist(-kRange, kRange - 1);
        std::uniform_int_distribution<int> mask_dist(0, 0xFF);

        PackedI64Words input{};
        for (std::size_t limb = 0; limb < input.size(); ++limb) {
            for (std::size_t lane = 0; lane < input[limb].size(); ++lane) {
                input[limb][lane] = dist(gen);
            }
        }

        const u8 lane_mask = static_cast<u8>(mask_dist(gen));
        PackedI64Words expected = input;
        const auto& mod = v2_scaled_words<Scale>();
        for (std::size_t limb = 0; limb < NWORDS; ++limb) {
            for (std::size_t lane = 0; lane < 8; ++lane) {
                if (((lane_mask >> lane) & 1u) != 0) {
                    expected[NWORDS + limb][lane] += static_cast<i64>(mod[limb]);
                }
            }
        }
        scalar_signed_carry_range(expected, NWORDS, (2 * NWORDS) - 1);
        const u8 expected_mask = high_negative_mask_scalar(expected);

        UFpPack8x1 got{};
        pack_from_i64_words(got, input);
        const u8 got_mask = lane_mask_to_u8(v2_sign_step<Scale>(got.data(), u8_to_lane_mask(lane_mask)));
        const PackedI64Words got_words = pack_to_i64_words(got);

        if (!packed_words_equal(got_words, expected) || got_mask != expected_mask) {
            std::fprintf(stderr, "v2_sign_step<%llu> validation failed\n", static_cast<unsigned long long>(Scale));
            std::fprintf(stderr, "lane_mask=%u expected_mask=%u got_mask=%u\n",
                static_cast<unsigned>(lane_mask),
                static_cast<unsigned>(expected_mask),
                static_cast<unsigned>(got_mask));
            dump_packed_words("expected", expected);
            dump_packed_words("got", got_words);
            std::abort();
        }
    }


    void random_redc_fp(Fp& v, std::mt19937_64& gen)
    {
        for (u64 i = 0; i < v.num_limb(); ++i) {
            v[i] = gen();
        }
        while (Fp::is_redc(v) == 0) {
            Fp::rshift(v, v, 1);
        }
    }

    BLST_FORCEINLINE void bench_consume_fp6_wide(const std::vector<UFp6Pack8x1>& data)
    {
        if (data.empty()) {
            return;
        }
        benchmark::DoNotOptimize(data.front().c[0][0][0]);
        benchmark::DoNotOptimize(data.back().c[2][1][NWORDS - 1]);
        benchmark::ClobberMemory();
    }

    Fp2 make_fp2(const Fp& re, const Fp& im)
    {
        Fp2 out{};
        out[0] = re;
        out[1] = im;
        return out;
    }

    Fp6 make_fp6(const Fp2& c0, const Fp2& c1, const Fp2& c2)
    {
        Fp6 out{};
        out[0] = c0;
        out[1] = c1;
        out[2] = c2;
        return out;
    }

    Fp6 adversarial_fp6_from_code(const std::array<Fp2, 9>& patterns, std::size_t code, const bool is_a)
    {
        std::array<std::size_t, 3> digits{};
        const std::size_t offset = is_a ? 0 : 3;
        for (std::size_t i = 0; i < 3; ++i) {
            digits[i] = code % patterns.size();
            code /= patterns.size();
        }
        return make_fp6(patterns[digits[0 + offset - offset]], patterns[digits[1 + offset - offset]], patterns[digits[2 + offset - offset]]);
    }

    void decode_adversarial_case(
        const std::array<Fp2, 9>& patterns, const std::size_t code, Fp6& a, Fp6& b)
    {
        std::size_t tmp = code;
        std::array<std::size_t, 6> digits{};
        for (std::size_t i = 0; i < digits.size(); ++i) {
            digits[i] = tmp % patterns.size();
            tmp /= patterns.size();
        }

        a = make_fp6(patterns[digits[0]], patterns[digits[1]], patterns[digits[2]]);
        b = make_fp6(patterns[digits[3]], patterns[digits[4]], patterns[digits[5]]);
    }

    BLST_FORCEINLINE void prepare_lufp_v2_fp6_raw(UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
    {
        mul_fp6_pack_lufpv2_wide<false>(out, a, b);
    }

    struct Fp6Pack8x1V2BenchData
    {
        explicit Fp6Pack8x1V2BenchData(const std::size_t n_packs)
            : n_packs_(n_packs)
            , a_scalar(n_packs_ * 8)
            , b_scalar(n_packs_ * 8)
            , a_pack(n_packs_)
            , b_pack(n_packs_)
            , out_wide(n_packs_)
        {
            std::mt19937_64 gen(0x8C4F7D2E19B3ULL + static_cast<u64>(n_packs_));

            for (std::size_t i = 0; i < n_packs_ * 8; ++i) {
                for (std::size_t c = 0; c < 3; ++c) {
                    for (std::size_t d = 0; d < 2; ++d) {
                        random_redc_fp(a_scalar[i][c][d], gen);
                        random_redc_fp(b_scalar[i][c][d], gen);
                    }
                }
            }

            for (std::size_t p = 0; p < n_packs_; ++p) {
                Fp6Pack8x1::from_aos(a_pack[p], a_scalar.data() + (8 * p));
                Fp6Pack8x1::from_aos(b_pack[p], b_scalar.data() + (8 * p));
            }
        }

        const std::size_t n_packs_;
        std::vector<Fp6> a_scalar;
        std::vector<Fp6> b_scalar;
        std::vector<Fp6Pack8x1> a_pack;
        std::vector<Fp6Pack8x1> b_pack;
        std::vector<UFp6Pack8x1> out_wide;
    };

    template <typename MulFn>
    void validate_lufp_v2_fp6_kernel(const char* label, MulFn&& mul_fn, bool validate_core_ops = true);

    template <typename MulFn>
    void validate_lufp_v2_fp6_kernel(const char* label, MulFn&& mul_fn, const bool /*validate_core_ops*/)
    {
        Fp6Pack8x1V2BenchData data(1);
        std::array<Fp6, 8> got_candidate{};
        std::array<Fp6, 8> ref{};

        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            Fp6Pack8x1 candidate_redc{};
            Fp6Pack8x1 ref_mul{};
            UFp6Pack8x1 candidate;

            mul_fn(candidate, data.a_pack[p], data.b_pack[p]);
            UFp6Pack8x1::redc(candidate_redc, candidate);
            Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

            Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
            Fp6Pack8x1::to_aos(ref.data(), ref_mul);

            for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
                if (!fp6_equal(got_candidate[lane], ref[lane])) {
                    std::fprintf(stderr, "Fp6 LUFp v2 %s validation failed at pack %zu lane %zu\n", label, p, lane);
                    dump_fp6("candidate", got_candidate[lane]);
                    dump_fp6("ref", ref[lane]);
                    std::abort();
                }
            }
        }
    }

    void validate_lufp_v2_fp6()
    {
        validate_lufp_v2_fp6_kernel(
            "wide",
            [](UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b) {
                mul_fp6_pack_lufpv2_wide(out, a, b);
            });
    }



    void renorm_negative_fold_fp6(UFp6Pack8x1& in_out)
    {
        lufp_v2_renorm_negative_fold<3u, 4u>(in_out.c[0][0]);
        lufp_v2_renorm_negative_fold<4u, 3u>(in_out.c[0][1]);
        lufp_v2_renorm_negative_fold<3u, 2u>(in_out.c[1][0]);
        lufp_v2_renorm_negative_fold<2u, 3u>(in_out.c[1][1]);
        lufp_v2_renorm_negative_fold<2u, 2u>(in_out.c[2][0]);
        lufp_v2_renorm_negative_fold<2u, 2u>(in_out.c[2][1]);
    }


    void BM_Fp6MulLufpWide(benchmark::State& state)
    {
        Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

        for (auto _ : state) {
            for (std::size_t p = 0; p < data.n_packs_; ++p) {
                mul_fp6_pack_lufpv2_wide(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
            }
            bench_consume_fp6_wide(data.out_wide);
        }

        state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
            static_cast<int64_t>(data.n_packs_));
    }

    void BM_Fp6MulLufpWideLazy(benchmark::State& state)
    {
        Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

        for (auto _ : state) {
            for (std::size_t p = 0; p < data.n_packs_; ++p) {
                mul_fp6_pack_lufpv2_wide<false>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
            }
            bench_consume_fp6_wide(data.out_wide);
        }

        state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
            static_cast<int64_t>(data.n_packs_));
    }


    void renorm_negative_fold_fp6_fused(UFp6Pack8x1& in_out)
    {
        lufp_v2_renorm_negative_fold_fused<3u, 4u>(in_out.c[0][0]);
        lufp_v2_renorm_negative_fold_fused<4u, 3u>(in_out.c[0][1]);
        lufp_v2_renorm_negative_fold_fused<3u, 2u>(in_out.c[1][0]);
        lufp_v2_renorm_negative_fold_fused<2u, 3u>(in_out.c[1][1]);
        lufp_v2_renorm_negative_fold_fused<2u, 2u>(in_out.c[2][0]);
        lufp_v2_renorm_negative_fold_fused<2u, 2u>(in_out.c[2][1]);
    }

    template <void (*RenormFn)(UFp6Pack8x1&)>
    void BM_Fp6RawRenormNegativeFoldFused(benchmark::State& state)
    {
        Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

        for (auto _ : state) {
            state.PauseTiming();
            for (std::size_t p = 0; p < data.n_packs_; ++p) {
                prepare_lufp_v2_fp6_raw(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
            }
            state.ResumeTiming();

            for (std::size_t p = 0; p < data.n_packs_; ++p) {
                RenormFn(data.out_wide[p]);
            }
            bench_consume_fp6_wide(data.out_wide);
        }

        state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
            static_cast<int64_t>(data.n_packs_));
    }

    BENCHMARK(BM_Fp6MulLufpWide)->Arg(256);
    BENCHMARK(BM_Fp6MulLufpWideLazy)->Arg(256);
    BENCHMARK_TEMPLATE(BM_Fp6RawRenormNegativeFoldFused, renorm_negative_fold_fp6_fused)->Arg(256);

}
