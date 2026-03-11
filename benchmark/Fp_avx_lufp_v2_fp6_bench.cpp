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

    v2_signed_carry_range(safe_state.data(), 0, (2 * NWORDS) - 1);
    const u8 neg_mask = lane_mask_to_u8(v2_signbit_mask(safe_state[(2 * NWORDS) - 1]));
    v2_high_add_mod_masked_u8(safe_state.data(), neg_mask, kV2P48X4.data());
    v2_signed_carry_range(safe_state.data(), NWORDS, (2 * NWORDS) - 1);
    const u8 neg_fix_mask = lane_mask_to_u8(v2_signbit_mask(safe_state[(2 * NWORDS) - 1]));
    v2_high_add_mod_masked_u8(safe_state.data(), neg_fix_mask, kV2P48X1.data());
    v2_signed_carry_range(safe_state.data(), NWORDS, (2 * NWORDS) - 1);
    const u8 ge_4p_mask = v2_high_ge_mod_mask_u8(safe_state.data(), kV2P48X4.data());
    const u8 ge_2p_mask = v2_high_ge_mod_mask_u8(safe_state.data(), kV2P48X2.data());
    const u8 ge_p_mask = v2_high_ge_p_mask_u8(safe_state.data());

    v2_high_sub_const(candidate_state.data(), kV2Bias3PPlus3Vec);
    const u8 sign0 = lane_mask_to_u8(v2_carry_range_with_sign(candidate_state.data(), 0, (2 * NWORDS) - 1));
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

void debug_search_half_repairs(const char* half_label, const UFpPack8x1& raw, const UFpPack8x1& want)
{
    UFpPack8x1 got{};

    got = raw;
    lufp_v2_fp6_mul_common_renorm_half_safe(got);
    report_half_repair_match(half_label, "safe", got, want);

    got = raw;
    lufp_v2_fp6_mul_common_renorm_half_fast(got);
    report_half_repair_match(half_label, "fast", got, want);

    got = raw;
    lufp_v2_fp6_mul_common_renorm_half_step_candidate(got);
    report_half_repair_match(half_label, "final_step", got, want);

    got = raw;
    lufp_v2_renorm_exact_add(got);
    report_half_repair_match(half_label, "exact_add", got, want);

    got = raw;
    lufp_v2_renorm_exact_add(got);
    lufp_v2_renorm_exact_add(got);
    report_half_repair_match(half_label, "exact_add_twice", got, want);

    got = raw;
    lufp_v2_renorm_relaxed2_step211(got);
    report_half_repair_match(half_label, "step211", got, want);

    got = raw;
    lufp_v2_fp6_mul_common_renorm_half_step_candidate(got);
    lufp_v2_renorm_exact_add(got);
    report_half_repair_match(half_label, "final_step_then_exact_add", got, want);

    got = raw;
    lufp_v2_renorm_relaxed2_step211(got);
    lufp_v2_renorm_exact_add(got);
    report_half_repair_match(half_label, "step211_then_exact_add", got, want);

    got = raw;
    lufp_v2_fp6_mul_common_renorm_half_step_candidate(got);
    lufp_v2_renorm_relaxed2_step211(got);
    report_half_repair_match(half_label, "final_step_then_step211", got, want);
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
            v2_signed_carry_range(got.data(), 0, (2 * NWORDS) - 1);
            v2_lane_mask_t sign_mask = v2_signbit_mask(got[(2 * NWORDS) - 1]);
            for (int i = 0; i < len; ++i) {
                const int digit = digits[static_cast<std::size_t>(i)];
                if (kScales[digit] == 1) {
                    sign_mask = v2_sign_step<1>(got.data(), sign_mask);
                } else if (kScales[digit] == 2) {
                    sign_mask = v2_sign_step<2>(got.data(), sign_mask);
                } else {
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

void maybe_renorm_exact_add(UFp2Pack8x1x1& in_out, const u32 mask, const u32 bit)
{
    if ((mask & (1u << bit)) != 0u) {
        lufp2_v2_renorm_exact_add(in_out);
    }
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

void debug_compare_step_candidate_manual_sequence(const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
    Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{};
    UFp2Pack8x1x1 w0{}, w1{}, w2{};

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

    auto compare_pair = [&](const char* base_label,
                            const UFp2Pack8x1x1& w,
                            const UFp2Pack8x1x1& s0,
                            const UFp2Pack8x1x1& s1,
                            const UFp2Pack8x1x1* extra_add)
    {
        UFp2Pack8x1x1 cand = w;
        UFp2Pack8x1x1 safe = w;
        UFp2Pack8x1x1 tmp_cand{};
        UFp2Pack8x1x1 tmp_safe{};

        lufp2_v2_sub(cand, cand, s0);
        lufp2_v2_sub(safe, safe, s0);
        UFp2Pack8x1x1 exact = cand;
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        lufp2_v2_renorm_exact_add_pair(exact);
        if (!fp2_pair_words_equal(cand, safe)) {
            std::fprintf(stderr, "%s exact-add match after first sub: %s\n",
                         base_label,
                         fp2_pair_words_equal(exact, safe) ? "yes" : "no");
        }
        assert_fp2_pair_equal_or_dump((std::string(base_label) + " after first sub").c_str(), cand, safe);

        lufp2_v2_sub(cand, cand, s1);
        lufp2_v2_sub(safe, safe, s1);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump((std::string(base_label) + " after second sub").c_str(), cand, safe);

        if (extra_add != nullptr) {
            lufp2_v2_mul_by_u_plus_1(tmp_cand, cand);
            lufp2_v2_mul_by_u_plus_1(tmp_safe, safe);
            lufp2_v2_renorm_fp6_mul_common_step_candidate(tmp_cand);
            lufp2_v2_renorm_fp6_mul_common_safe(tmp_safe);
            assert_fp2_pair_equal_or_dump((std::string(base_label) + " after u+1").c_str(), tmp_cand, tmp_safe);

            lufp2_v2_add(cand, tmp_cand, *extra_add);
            lufp2_v2_add(safe, tmp_safe, *extra_add);
            lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
            lufp2_v2_renorm_fp6_mul_common_safe(safe);
            assert_fp2_pair_equal_or_dump((std::string(base_label) + " after final add").c_str(), cand, safe);
        }
    };

    compare_pair("c0", w0, t1, t2, &t0);

    {
        UFp2Pack8x1x1 cand = w1;
        UFp2Pack8x1x1 safe = w1;
        UFp2Pack8x1x1 u1_cand{};
        UFp2Pack8x1x1 u1_safe{};

        lufp2_v2_sub(cand, cand, t0);
        lufp2_v2_sub(safe, safe, t0);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump("c1 after first sub", cand, safe);

        lufp2_v2_sub(cand, cand, t1);
        lufp2_v2_sub(safe, safe, t1);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump("c1 after second sub", cand, safe);

        lufp2_v2_mul_by_u_plus_1(u1_cand, t2);
        lufp2_v2_mul_by_u_plus_1(u1_safe, t2);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(u1_cand);
        lufp2_v2_renorm_fp6_mul_common_safe(u1_safe);
        assert_fp2_pair_equal_or_dump("c1 after u+1", u1_cand, u1_safe);

        lufp2_v2_add(cand, cand, u1_cand);
        lufp2_v2_add(safe, safe, u1_safe);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump("c1 after final add", cand, safe);
    }

    {
        UFp2Pack8x1x1 cand = w2;
        UFp2Pack8x1x1 safe = w2;

        lufp2_v2_sub(cand, cand, t0);
        lufp2_v2_sub(safe, safe, t0);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump("c2 after first sub", cand, safe);

        lufp2_v2_sub(cand, cand, t2);
        lufp2_v2_sub(safe, safe, t2);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump("c2 after second sub", cand, safe);

        lufp2_v2_add(cand, cand, t1);
        lufp2_v2_add(safe, safe, t1);
        lufp2_v2_renorm_fp6_mul_common_step_candidate(cand);
        lufp2_v2_renorm_fp6_mul_common_safe(safe);
        assert_fp2_pair_equal_or_dump("c2 after final add", cand, safe);
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

void validate_lufp_v2_core_ops();

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

template <typename MulFn>
void validate_lufp_v2_fp6_adversarial(const char* label, MulFn&& mul_fn)
{
    validate_lufp_v2_core_ops();

    const Fp zero = Fp::zero();
    const Fp one = Fp::one();
    Fp minus_one{};
    Fp::neg(minus_one, one);

    const std::array<Fp2, 9> patterns = {
        make_fp2(zero, zero),
        make_fp2(one, zero),
        make_fp2(zero, one),
        make_fp2(one, one),
        make_fp2(minus_one, zero),
        make_fp2(zero, minus_one),
        make_fp2(minus_one, minus_one),
        make_fp2(one, minus_one),
        make_fp2(minus_one, one),
    };

    constexpr std::size_t kBatch = 8;
    constexpr std::size_t kTotalCases = 531441; // 9^6

    std::array<Fp6, kBatch> a_scalar{};
    std::array<Fp6, kBatch> b_scalar{};
    std::array<Fp6, kBatch> got_candidate{};
    std::array<Fp6, kBatch> got_safe{};
    std::array<Fp6, kBatch> ref{};

    for (std::size_t base = 0; base < kTotalCases; base += kBatch) {
        for (std::size_t lane = 0; lane < kBatch; ++lane) {
            const std::size_t case_idx = (base + lane < kTotalCases) ? (base + lane) : (kTotalCases - 1);
            decode_adversarial_case(patterns, case_idx, a_scalar[lane], b_scalar[lane]);
        }

        Fp6Pack8x1 a_pack{};
        Fp6Pack8x1 b_pack{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};

        Fp6Pack8x1::from_aos(a_pack, a_scalar.data());
        Fp6Pack8x1::from_aos(b_pack, b_scalar.data());

        mul_fn(candidate, a_pack, b_pack);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, a_pack, b_pack);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, a_pack, b_pack);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < kBatch && (base + lane) < kTotalCases; ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr,
                             "%s adversarial validation failed at case %zu lane %zu\n",
                             label,
                             base + lane,
                             lane);
                dump_fp6("a", a_scalar[lane]);
                dump_fp6("b", b_scalar[lane]);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }

    std::fprintf(stderr, "%s adversarial validation checked %zu structured cases\n", label, kTotalCases);
}

void validate_lufp_v2_fp6_step_candidate_sparse5b_adversarial()
{
    validate_lufp_v2_fp6_adversarial(
        "Sparse5B",
        [](UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b) {
            mul_fp6_pack_lufpv2_wide_step_candidate_sparse<0x3C2u>(out, a, b);
        });
}

void validate_lufp_v2_fp6_delayed_block_final_adversarial()
{
    validate_lufp_v2_fp6_adversarial(
        "DelayedBlockFinal",
        [](UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b) {
            mul_fp6_pack_lufpv2_wide_delayed_block_final(out, a, b);
        });
}

//struct LUFp6Pack8x1Raw
//{
//    LUFp2Pack8x1x1 c[3];
//};

BLST_FORCEINLINE void prepare_lufp_v2_fp6_raw(UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
    Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{};
    UFp2Pack8x1x1 w0{}, w1{}, w2{};

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

    mul_fp6_pack_lufpv2_combine_only_raw(
        out.c[0],
        out.c[1],
        out.c[2],
        t0,
        t1,
        t2,
        w0,
        w1,
        w2);
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
void validate_lufp_v2_fp6_kernel(const char* label, MulFn&& mul_fn);

void validate_lufp_v2_core_ops()
{
    std::mt19937_64 gen(0xB17F6A2D53ULL);
    constexpr i64 kRange = static_cast<i64>(1ULL << 47);
    std::uniform_int_distribution<i64> dist(-kRange, kRange - 1);

    {
        PackedI64Words a_words{}, b_words{}, c_words{}, expected{};
        for (std::size_t limb = 0; limb < expected.size(); ++limb) {
            for (std::size_t lane = 0; lane < expected[limb].size(); ++lane) {
                a_words[limb][lane] = dist(gen);
                b_words[limb][lane] = dist(gen);
                c_words[limb][lane] = dist(gen);
                expected[limb][lane] = a_words[limb][lane] - b_words[limb][lane] - c_words[limb][lane];
            }
        }

        UFpPack8x1 a{}, b{}, c{}, got{};
        pack_from_i64_words(a, a_words);
        pack_from_i64_words(b, b_words);
        pack_from_i64_words(c, c_words);
        lufp_v2_sub2(got, a, b, c);

        const PackedI64Words got_words = pack_to_i64_words(got);
        if (!packed_words_equal(got_words, expected)) {
            std::fprintf(stderr, "lufp_v2_sub2 validation failed\n");
            dump_packed_words("expected", expected);
            dump_packed_words("got", got_words);
            std::abort();
        }
    }

    {
        PackedI64Words input{};
        for (std::size_t limb = 0; limb < input.size(); ++limb) {
            for (std::size_t lane = 0; lane < input[limb].size(); ++lane) {
                input[limb][lane] = dist(gen) * 3;
            }
        }

        PackedI64Words expected = input;
        scalar_signed_carry_range(expected, 0, (2 * NWORDS) - 1);

        UFpPack8x1 got{};
        pack_from_i64_words(got, input);
        v2_signed_carry_range(got.data(), 0, (2 * NWORDS) - 1);
        const PackedI64Words got_words = pack_to_i64_words(got);

        if (!packed_words_equal(got_words, expected)) {
            std::fprintf(stderr, "v2_signed_carry_range validation failed\n");
            dump_packed_words("expected", expected);
            dump_packed_words("got", got_words);
            std::abort();
        }
    }

    validate_v2_sign_step_case<1>(gen);
    validate_v2_sign_step_case<2>(gen);
    validate_v2_sign_step_case<4>(gen);

    {
        Fp6Pack8x1V2BenchData data(2);
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            UFp6Pack8x1 raw{};
            prepare_lufp_v2_fp6_raw(raw, data.a_pack[p], data.b_pack[p]);

            for (std::size_t c = 0; c < 3; ++c) {
                for (std::size_t d = 0; d < 2; ++d) {
                    UFpPack8x1 safe = raw.c[c][d];
                    lufp_v2_fp6_mul_common_renorm_half_safe(safe);

                    LUFpPack8x1 legacy{};
                    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
                        legacy[limb] = raw.c[c][d][limb];
                    }

                    UFpPack8x1 ref{};
                    LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg4Pos4>(ref, legacy);

                    const PackedI64Words safe_words = pack_to_i64_words(safe);
                    const PackedI64Words ref_words = pack_to_i64_words(ref);
                    if (!packed_words_equal(safe_words, ref_words)) {
                        std::fprintf(stderr, "lufp_v2_fp6_mul_common_renorm_half_safe validation failed at pack %zu c%zu.%s\n",
                                     p,
                                     c,
                                     d == 0 ? "re" : "im");
                        dump_packed_words("safe", safe_words);
                        dump_packed_words("ref", ref_words);
                        std::abort();
                    }
                }
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

void validate_lufp_v2_fp6_step_candidate()
{
    validate_lufp_v2_core_ops();
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got_candidate{};
    std::array<Fp6, 8> got_safe{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};

        mul_fp6_pack_lufpv2_wide_step_candidate(candidate, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 step-candidate validation failed at pack %zu lane %zu\n", p, lane);
                debug_compare_step_candidate_manual_sequence(data.a_pack[p], data.b_pack[p]);
                UFp6Pack8x1 raw{};
                prepare_lufp_v2_fp6_raw(raw, data.a_pack[p], data.b_pack[p]);
                dump_step_candidate_half_debug("c0.re", raw.c[0][0]);
                dump_step_candidate_half_debug("c0.im", raw.c[0][1]);
                dump_step_candidate_half_debug("c1.re", raw.c[1][0]);
                dump_step_candidate_half_debug("c1.im", raw.c[1][1]);
                dump_step_candidate_half_debug("c2.re", raw.c[2][0]);
                dump_step_candidate_half_debug("c2.im", raw.c[2][1]);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }
}

template <typename MulFn>
void validate_lufp_v2_fp6_kernel(const char* label, MulFn&& mul_fn)
{
    validate_lufp_v2_core_ops();
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got_candidate{};
    std::array<Fp6, 8> got_safe{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};

        mul_fn(candidate, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 %s validation failed at pack %zu lane %zu\n", label, p, lane);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }
}

void validate_lufp_v2_fp6_exact_add_every_op()
{
    validate_lufp_v2_fp6_kernel(
        "exact-add-every-op",
        [](UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b) {
            mul_fp6_pack_lufpv2_wide_exact_add_every_op(out, a, b);
        });
}

template <u32 Mask>
void validate_lufp_v2_fp6_step_candidate_sparse(const char* label)
{
    validate_lufp_v2_core_ops();
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got_candidate{};
    std::array<Fp6, 8> got_safe{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};

        mul_fp6_pack_lufpv2_wide_step_candidate_sparse<Mask>(candidate, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 %s validation failed at pack %zu lane %zu\n", label, p, lane);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }
}

template <u32 ReMask, u32 ImMask>
void validate_lufp_v2_fp6_step_candidate_half_sparse(const char* label)
{
    validate_lufp_v2_core_ops();
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got_candidate{};
    std::array<Fp6, 8> got_safe{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};

        mul_fp6_pack_lufpv2_wide_step_candidate_half_sparse<ReMask, ImMask>(candidate, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 %s validation failed at pack %zu lane %zu\n", label, p, lane);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }
}

void validate_lufp_v2_fp6_step_candidate_relaxed2()
{
    validate_lufp_v2_core_ops();
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got_candidate{};
    std::array<Fp6, 8> got_safe{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};

        mul_fp6_pack_lufpv2_wide_step_candidate_relaxed2(candidate, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 step-candidate-relaxed2 validation failed at pack %zu lane %zu\n", p, lane);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }
}

void validate_lufp_v2_fp6_step_candidate_hybrid()
{
    validate_lufp_v2_core_ops();
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got_candidate{};
    std::array<Fp6, 8> got_safe{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 candidate_redc{};
        Fp6Pack8x1 safe_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 candidate{};
        UFp6Pack8x1 safe{};

        mul_fp6_pack_lufpv2_wide_step_candidate_hybrid(candidate, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(candidate_redc, candidate);

        mul_fp6_pack_lufpv2_wide_safe(safe, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(safe_redc, safe);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got_candidate.data(), candidate_redc);
        Fp6Pack8x1::to_aos(got_safe.data(), safe_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);

        for (std::size_t lane = 0; lane < got_candidate.size(); ++lane) {
            if (!fp6_equal(got_candidate[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 step-candidate-hybrid validation failed at pack %zu lane %zu\n", p, lane);
                dump_fp6("candidate", got_candidate[lane]);
                dump_fp6("safe", got_safe[lane]);
                dump_fp6("ref", ref[lane]);
                std::abort();
            }
        }
    }
}

enum class LUFpFp6MulStage
{
    PreAdd,
    PreAddMul,
    PreAddMulRawCombine,
    Full,
};

template <LUFpFp6MulStage Stage>
BLST_FORCEINLINE void lufp_fp6_mul_accounting_v2(Fp6Pack8x1V2BenchData& data)
{
    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
        Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
        UFp2Pack8x1x1 t0{}, t1{}, t2{};
        UFp2Pack8x1x1 w0{}, w1{}, w2{};

        Fp2Pack8x1x1::add(aa0, data.a_pack[p].c[1], data.a_pack[p].c[2]);
        Fp2Pack8x1x1::add(bb0, data.b_pack[p].c[1], data.b_pack[p].c[2]);
        Fp2Pack8x1x1::add(aa1, data.a_pack[p].c[0], data.a_pack[p].c[1]);
        Fp2Pack8x1x1::add(bb1, data.b_pack[p].c[0], data.b_pack[p].c[1]);
        Fp2Pack8x1x1::add(aa2, data.a_pack[p].c[0], data.a_pack[p].c[2]);
        Fp2Pack8x1x1::add(bb2, data.b_pack[p].c[0], data.b_pack[p].c[2]);

        if constexpr (Stage >= LUFpFp6MulStage::PreAddMul) {
            UFp2Pack8x1x1::mul(t0, data.a_pack[p].c[0], data.b_pack[p].c[0]);
            UFp2Pack8x1x1::mul(t1, data.a_pack[p].c[1], data.b_pack[p].c[1]);
            UFp2Pack8x1x1::mul(t2, data.a_pack[p].c[2], data.b_pack[p].c[2]);
            UFp2Pack8x1x1::mul(w0, aa0, bb0);
            UFp2Pack8x1x1::mul(w1, aa1, bb1);
            UFp2Pack8x1x1::mul(w2, aa2, bb2);

            if constexpr (Stage >= LUFpFp6MulStage::PreAddMulRawCombine) {
                if constexpr (Stage >= LUFpFp6MulStage::Full) {
                    mul_fp6_pack_lufpv2_combine_only_raw(
                        data.out_wide[p].c[0],
                        data.out_wide[p].c[1],
                        data.out_wide[p].c[2],
                        t0,
                        t1,
                        t2,
                        w0,
                        w1,
                        w2);
                    lufp2_v2_renorm_fp6_mul_common_safe(data.out_wide[p].c[0]);
                    lufp2_v2_renorm_fp6_mul_common_safe(data.out_wide[p].c[1]);
                    lufp2_v2_renorm_fp6_mul_common_safe(data.out_wide[p].c[2]);
                } else {
                    UFp2Pack8x1x1 c0{}, c1{}, c2{};
                    mul_fp6_pack_lufpv2_combine_only_raw(
                        c0,
                        c1,
                        c2,
                        t0,
                        t1,
                        t2,
                        w0,
                        w1,
                        w2);
                    benchmark::DoNotOptimize(c0[0][0]);
                    benchmark::DoNotOptimize(c2[1][NWORDS - 1]);
                }
            } else {
                benchmark::DoNotOptimize(t0[0][0]);
                benchmark::DoNotOptimize(w2[1][NWORDS - 1]);
            }
        } else {
            benchmark::DoNotOptimize(aa0[0][0]);
            benchmark::DoNotOptimize(bb2[1][NWORDS - 1]);
        }
    }
}

template <LUFpFp6MulStage Stage>
void BM_Fp6MulLufpAccounting(benchmark::State& state)
{
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        lufp_fp6_mul_accounting_v2<Stage>(data);
        if constexpr (Stage >= LUFpFp6MulStage::Full) {
            bench_consume_fp6_wide(data.out_wide);
        } else {
            benchmark::ClobberMemory();
        }
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

enum class LUFpV2Fp6MulCommonRenormStage
{
    Bias,
    BiasCarry,
    BiasCarrySign4,
    BiasCarrySign42,
    BiasCarrySign421,
    BiasCarrySign4211,
    Full,
};

template <LUFpV2Fp6MulCommonRenormStage Stage>
BLST_FORCEINLINE void lufp_v2_fp6_mul_common_renorm_accounting(UFpPack8x1& in_out)
{
    v2_high_sub_const(in_out.data(), kV2Bias3PPlus3Vec);

    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarry) {
        u8 sign_mask = v2_carry_range_with_sign(in_out.data(), 0, (2 * NWORDS) - 1);

        if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4) {
            sign_mask = v2_sign_step<4>(in_out.data(), sign_mask);

            if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign42) {
                sign_mask = v2_sign_step<2>(in_out.data(), sign_mask);

                if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign421) {
                    sign_mask = v2_sign_step<1>(in_out.data(), sign_mask);

                    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4211) {
                        benchmark::DoNotOptimize(v2_sign_step<1>(in_out.data(), sign_mask));
                    }
                }
            }
        }
    }

    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
        benchmark::DoNotOptimize(in_out[0]);
        benchmark::DoNotOptimize(in_out[(2 * NWORDS) - 1]);
    } else {
        benchmark::DoNotOptimize(in_out[0]);
        benchmark::DoNotOptimize(in_out[(2 * NWORDS) - 1]);
    }
}

template <LUFpV2Fp6MulCommonRenormStage Stage>
void BM_LUFpV2Fp6MulCommonRenormAccounting(benchmark::State& state)
{
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        state.PauseTiming();
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            prepare_lufp_v2_fp6_raw(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        state.ResumeTiming();

        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[0][0]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[0][1]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[1][0]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[1][1]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[2][0]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[2][1]);
        }

        benchmark::ClobberMemory();
        if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
            bench_consume_fp6_wide(data.out_wide);
        }
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(6));
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

void BM_Fp6MulLufpWideSafe(benchmark::State& state)
{
    validate_lufp_v2_fp6();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_safe(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidate(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideExactAddEveryOp(benchmark::State& state)
{
    validate_lufp_v2_fp6_exact_add_every_op();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_exact_add_every_op(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateSparse5A(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_sparse<0x362u>("step-candidate-sparse5a");
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_sparse<0x362u>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateSparse5B(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_sparse<0x3C2u>("step-candidate-sparse5b");
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_sparse<0x3C2u>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateSparse5C(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_sparse<0x562u>("step-candidate-sparse5c");
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_sparse<0x562u>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateSparse5D(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_sparse<0x5C2u>("step-candidate-sparse5d");
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_sparse<0x5C2u>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateHalfSparse8A(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_half_sparse<0x5C2u, 0x482u>("step-candidate-half-sparse8a");
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_half_sparse<0x5C2u, 0x482u>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateHalfSparse8B(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_half_sparse<0x5C2u, 0x182u>("step-candidate-half-sparse8b");
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_half_sparse<0x5C2u, 0x182u>(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateRelaxed2(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_relaxed2();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_relaxed2(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateHybrid(benchmark::State& state)
{
    validate_lufp_v2_fp6_step_candidate_hybrid();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide_step_candidate_hybrid(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWideStepCandidateSparse5BAdversarialSearch(benchmark::State& state)
{
    static bool validated = false;
    if (!validated) {
        validate_lufp_v2_fp6_step_candidate_sparse5b_adversarial();
        validated = true;
    }

    for (auto _ : state) {
        benchmark::DoNotOptimize(validated);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

void BM_Fp6MulLufpWideDelayedBlockFinalAdversarialSearch(benchmark::State& state)
{
    static bool validated = false;
    if (!validated) {
        validate_lufp_v2_fp6_delayed_block_final_adversarial();
        validated = true;
    }

    for (auto _ : state) {
        benchmark::DoNotOptimize(validated);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

BENCHMARK(BM_Fp6MulLufpWide)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideSafe)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidate)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideExactAddEveryOp)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateSparse5A)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateSparse5B)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateSparse5C)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateSparse5D)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateHalfSparse8A)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateHalfSparse8B)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateRelaxed2)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateHybrid)->Arg(256);
BENCHMARK(BM_Fp6MulLufpWideStepCandidateSparse5BAdversarialSearch);
BENCHMARK(BM_Fp6MulLufpWideDelayedBlockFinalAdversarialSearch);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMul)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMulRawCombine)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::Full)->Arg(256);

BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::Bias)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign42)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign421)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4211)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::Full)->Arg(256);
}
