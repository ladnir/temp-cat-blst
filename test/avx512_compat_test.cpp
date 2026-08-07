#include "blst/common/defines.h"

#if !defined(BLST_AVX512) || !defined(BLST_AVX512_PORTABLE)
#error "avx512_compat_test requires BLST_AVX512 and BLST_AVX512_PORTABLE"
#endif

#include "blst/backend/avx512/intrin_compat.h"

#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <random>
#include <type_traits>

namespace blst::test
{
#if defined(__GNUC__) || defined(__clang__)
    static_assert(__is_same(vv512i, ::__m512i));
    static_assert(__is_same(vv256i, ::__m256i));
    static_assert(!__is_same(cc512i, ::__m512i));
#else
    static_assert(std::is_same_v<vv512i, ::__m512i>);
    static_assert(std::is_same_v<vv256i, ::__m256i>);
    static_assert(!std::is_same_v<cc512i, ::__m512i>);
#endif
    static_assert(sizeof(cc512i) == 64);
    static_assert(sizeof(cc256i) == 32);
    static_assert(sizeof(cc128i) == 16);
    static_assert(std::is_standard_layout_v<cc512i>);
    static_assert(std::is_trivially_copyable_v<cc512i>);

    class Avx512CompatTest : public ::testing::Test
    {
    protected:
        std::mt19937_64 rng{0x434154424c535437ULL};

        auto random_cc512() -> cc512i
        {
            cc512i out{};
            for (u64& lane : out.lane) {
                lane = rng();
            }
            return out;
        }

        static void expect_equal(const ::__m512i expected, const cc512i actual)
        {
            const cc512i expected_cc = _cc512_from_mm512(expected);
            for (u64 i = 0; i < 8; ++i) {
                EXPECT_EQ(actual.lane[i], expected_cc.lane[i]) << "lane " << i;
            }
        }

        static void expect_equal(const ::__m256i expected, const cc256i actual)
        {
            const cc256i expected_cc = _cc256_from_mm256(expected);
            for (u64 i = 0; i < 4; ++i) {
                EXPECT_EQ(actual.lane[i], expected_cc.lane[i]) << "lane " << i;
            }
        }

        static auto mm(const cc512i value) -> ::__m512i
        {
            return _mm512_from_cc512(value);
        }

        static void expect_equal(const ::__m128i expected, const cc128i actual)
        {
            const cc128i expected_cc = _cc128_from_mm128(expected);
            EXPECT_EQ(actual.lane[0], expected_cc.lane[0]);
            EXPECT_EQ(actual.lane[1], expected_cc.lane[1]);
        }

        template <int bits>
        static void check_immediate_shifts(const cc512i a)
        {
            expect_equal(_mm512_slli_epi64(mm(a), bits), _cc512_slli_epi64(a, bits));
            expect_equal(_mm512_srli_epi64(mm(a), bits), _cc512_srli_epi64(a, bits));
            expect_equal(_mm512_srai_epi64(mm(a), bits), _cc512_srai_epi64(a, bits));
        }

        template <int imm8>
        static void check_immediate_permutations(const cc512i a, u64 mask8, u64 mask16)
        {
            expect_equal(
                _mm512_shuffle_epi32(mm(a), static_cast<_MM_PERM_ENUM>(imm8)),
                _cc512_shuffle_epi32<imm8>(a));
            expect_equal(
                _mm512_maskz_shuffle_epi32(static_cast<__mmask16>(mask16), mm(a), static_cast<_MM_PERM_ENUM>(imm8)),
                _cc512_maskz_shuffle_epi32<imm8>(mask16, a));
            expect_equal(
                _mm512_permutex_epi64(mm(a), imm8),
                _cc512_permutex_epi64<imm8>(a));
            expect_equal(
                _mm512_maskz_permutex_epi64(static_cast<__mmask8>(mask8), mm(a), imm8),
                _cc512_maskz_permutex_epi64<imm8>(mask8, a));
        }

        static void check_madd52_case(const cc512i acc, const cc512i a, const cc512i b)
        {
            expect_equal(
                _mm512_madd52lo_epu64(mm(acc), mm(a), mm(b)),
                _cc512_madd52lo_epu64(acc, a, b));
            expect_equal(
                _mm512_madd52hi_epu64(mm(acc), mm(a), mm(b)),
                _cc512_madd52hi_epu64(acc, a, b));
        }
    };

    TEST_F(Avx512CompatTest, Construction)
    {
        expect_equal(_mm512_setzero_si512(), _cc512_setzero_si512());

        constexpr std::array<u64, 5> values = {
            0, 1, u64{1} << 51, (u64{1} << 52) - 1, std::numeric_limits<u64>::max()
        };
        for (const u64 value : values) {
            expect_equal(
                _mm512_set1_epi64(static_cast<long long>(value)),
                _cc512_set1_epi64(static_cast<long long>(value)));
        }

        const std::array<u64, 8> lanes = {
            0, 1, 2, 3, std::numeric_limits<u64>::max(), u64{1} << 63, u64{1} << 52, 0x123456789abcdef0ULL
        };
        expect_equal(
            _mm512_set_epi64(
                static_cast<long long>(lanes[7]), static_cast<long long>(lanes[6]),
                static_cast<long long>(lanes[5]), static_cast<long long>(lanes[4]),
                static_cast<long long>(lanes[3]), static_cast<long long>(lanes[2]),
                static_cast<long long>(lanes[1]), static_cast<long long>(lanes[0])),
            _cc512_set_epi64(
                static_cast<long long>(lanes[7]), static_cast<long long>(lanes[6]),
                static_cast<long long>(lanes[5]), static_cast<long long>(lanes[4]),
                static_cast<long long>(lanes[3]), static_cast<long long>(lanes[2]),
                static_cast<long long>(lanes[1]), static_cast<long long>(lanes[0])));
    }

    TEST_F(Avx512CompatTest, ArithmeticLogicAndComparisons)
    {
        for (int iteration = 0; iteration < 128; ++iteration) {
            const cc512i a = random_cc512();
            const cc512i b = random_cc512();

            expect_equal(_mm512_add_epi64(mm(a), mm(b)), _cc512_add_epi64(a, b));
            expect_equal(_mm512_sub_epi64(mm(a), mm(b)), _cc512_sub_epi64(a, b));
            expect_equal(_mm512_and_si512(mm(a), mm(b)), _cc512_and_si512(a, b));
            expect_equal(_mm512_or_si512(mm(a), mm(b)), _cc512_or_si512(a, b));
            expect_equal(_mm512_xor_si512(mm(a), mm(b)), _cc512_xor_si512(a, b));
            expect_equal(_mm512_min_epu64(mm(a), mm(b)), _cc512_min_epu64(a, b));
            EXPECT_EQ(
                static_cast<u64>(_mm512_cmpgt_epi64_mask(mm(a), mm(b))),
                _cc512_cmpgt_epi64_mask(a, b));
            EXPECT_EQ(
                static_cast<u64>(_mm512_cmpge_epi64_mask(mm(a), mm(b))),
                _cc512_cmpge_epi64_mask(a, b));
        }
    }

    TEST_F(Avx512CompatTest, Shifts)
    {
        for (int iteration = 0; iteration < 32; ++iteration) {
            const cc512i a = random_cc512();
            check_immediate_shifts<0>(a);
            check_immediate_shifts<1>(a);
            check_immediate_shifts<26>(a);
            check_immediate_shifts<51>(a);
            check_immediate_shifts<63>(a);

            cc512i counts = random_cc512();
            constexpr std::array<u64, 8> boundary_counts = {0, 1, 26, 51, 63, 64, 65, std::numeric_limits<u64>::max()};
            if (iteration == 0) {
                counts.lane[0] = boundary_counts[0];
                counts.lane[1] = boundary_counts[1];
                counts.lane[2] = boundary_counts[2];
                counts.lane[3] = boundary_counts[3];
                counts.lane[4] = boundary_counts[4];
                counts.lane[5] = boundary_counts[5];
                counts.lane[6] = boundary_counts[6];
                counts.lane[7] = boundary_counts[7];
            }
            expect_equal(_mm512_sllv_epi64(mm(a), mm(counts)), _cc512_sllv_epi64(a, counts));
        }
    }

    TEST_F(Avx512CompatTest, LoadsAndStores)
    {
        alignas(64) std::array<u64, 8> aligned{};
        for (u64& lane : aligned) {
            lane = rng();
        }

        expect_equal(_mm512_load_si512(aligned.data()), _cc512_load_si512(aligned.data()));
        expect_equal(_mm512_loadu_si512(aligned.data()), _cc512_loadu_si512(aligned.data()));

        const cc512i value = random_cc512();
        alignas(64) std::array<u64, 8> cc_aligned{};
        alignas(64) std::array<u64, 8> mm_aligned{};
        _cc512_store_si512(cc_aligned.data(), value);
        _mm512_store_si512(mm_aligned.data(), mm(value));
        EXPECT_EQ(cc_aligned, mm_aligned);

        std::array<u64, 9> cc_unaligned{};
        std::array<u64, 9> mm_unaligned{};
        _cc512_storeu_si512(cc_unaligned.data() + 1, value);
        _mm512_storeu_si512(mm_unaligned.data() + 1, mm(value));
        EXPECT_EQ(cc_unaligned, mm_unaligned);
    }

    TEST_F(Avx512CompatTest, Masks)
    {
        constexpr std::array<u64, 5> masks = {0x00, 0x01, 0x55, 0xaa, 0xff};
        for (const u64 mask : masks) {
            const cc512i src = random_cc512();
            const cc512i a = random_cc512();
            const cc512i b = random_cc512();
            expect_equal(
                _mm512_mask_add_epi64(mm(src), static_cast<__mmask8>(mask), mm(a), mm(b)),
                _cc512_mask_add_epi64(src, mask, a, b));
            expect_equal(
                _mm512_mask_sub_epi64(mm(src), static_cast<__mmask8>(mask), mm(a), mm(b)),
                _cc512_mask_sub_epi64(src, mask, a, b));
            expect_equal(
                _mm512_mask_blend_epi64(static_cast<__mmask8>(mask), mm(a), mm(b)),
                _cc512_mask_blend_epi64(mask, a, b));
        }
    }

    TEST_F(Avx512CompatTest, Permutations)
    {
        for (int iteration = 0; iteration < 64; ++iteration) {
            const cc512i a = random_cc512();
            const cc512i indices = random_cc512();
            expect_equal(
                _mm512_permutexvar_epi64(mm(indices), mm(a)),
                _cc512_permutexvar_epi64(indices, a));
        }

        const cc512i a = random_cc512();
        check_immediate_permutations<0x00>(a, 0x00, 0x0000);
        check_immediate_permutations<0x1b>(a, 0x55, 0x5a5a);
        check_immediate_permutations<0x4e>(a, 0xaa, 0xa5a5);
        check_immediate_permutations<0xff>(a, 0xff, 0xffff);
    }

    TEST_F(Avx512CompatTest, Extraction)
    {
        const cc512i a = random_cc512();
        const cc256i low = _cc512_extracti64x4_epi64<0>(a);
        const cc256i high = _cc512_extracti64x4_epi64<1>(a);
        const ::__m256i mm_low = _mm512_extracti64x4_epi64(mm(a), 0);
        const ::__m256i mm_high = _mm512_extracti64x4_epi64(mm(a), 1);

        expect_equal(mm_low, low);
        expect_equal(mm_high, high);
        const ::__m256i mm_low_from_cc = _mm256_from_cc256(low);
        const ::__m256i mm_high_from_cc = _mm256_from_cc256(high);
        EXPECT_EQ(_cc256_extract_epi64<0>(low), static_cast<u64>(_mm256_extract_epi64(mm_low_from_cc, 0)));
        EXPECT_EQ(_cc256_extract_epi64<1>(low), static_cast<u64>(_mm256_extract_epi64(mm_low_from_cc, 1)));
        EXPECT_EQ(_cc256_extract_epi64<2>(high), static_cast<u64>(_mm256_extract_epi64(mm_high_from_cc, 2)));
        EXPECT_EQ(_cc256_extract_epi64<3>(high), static_cast<u64>(_mm256_extract_epi64(mm_high_from_cc, 3)));

        const cc128i lowest = _cc512_castsi512_si128(a);
        const ::__m128i mm_lowest = _mm512_castsi512_si128(mm(a));
        expect_equal(mm_lowest, lowest);
        EXPECT_EQ(_cc128_extract_epi64<0>(lowest), static_cast<u64>(_mm_cvtsi128_si64(mm_lowest)));
        EXPECT_EQ(_cc128_extract_epi64<1>(lowest), static_cast<u64>(_mm_extract_epi64(mm_lowest, 1)));
    }

    TEST_F(Avx512CompatTest, Low128HorizontalAndBroadcast)
    {
        constexpr std::array<int, 7> shifts = {0, 1, 26, 51, 52, 63, 64};
        for (int iteration = 0; iteration < 64; ++iteration) {
            const cc512i a = random_cc512();
            const cc512i b = random_cc512();
            const cc128i a128 = _cc512_castsi512_si128(a);
            const cc128i b128 = _cc512_castsi512_si128(b);
            const ::__m128i mm_a128 = _mm512_castsi512_si128(mm(a));
            const ::__m128i mm_b128 = _mm512_castsi512_si128(mm(b));

            expect_equal(_mm_add_epi64(mm_a128, mm_b128), _cc128_add_epi64(a128, b128));
            for (const int bits : shifts) {
                switch (bits) {
                case 0: expect_equal(_mm_srli_epi64(mm_a128, 0), _cc128_srli_epi64(a128, bits)); break;
                case 1: expect_equal(_mm_srli_epi64(mm_a128, 1), _cc128_srli_epi64(a128, bits)); break;
                case 26: expect_equal(_mm_srli_epi64(mm_a128, 26), _cc128_srli_epi64(a128, bits)); break;
                case 51: expect_equal(_mm_srli_epi64(mm_a128, 51), _cc128_srli_epi64(a128, bits)); break;
                case 52: expect_equal(_mm_srli_epi64(mm_a128, 52), _cc128_srli_epi64(a128, bits)); break;
                case 63: expect_equal(_mm_srli_epi64(mm_a128, 63), _cc128_srli_epi64(a128, bits)); break;
                case 64: expect_equal(_mm_srli_epi64(mm_a128, 64), _cc128_srli_epi64(a128, bits)); break;
                }
            }
            expect_equal(_mm512_broadcastq_epi64(mm_a128), _cc512_broadcastq_epi64(a128));
        }
    }

    TEST_F(Avx512CompatTest, Madd52BoundaryAndMaskedInputs)
    {
        constexpr u64 mask52 = (u64{1} << 52) - 1;
        const cc512i acc = {{
            0, 1, std::numeric_limits<u64>::max(), std::numeric_limits<u64>::max() - mask52,
            u64{1} << 63, mask52, mask52 + 1, 0xaaaaaaaaaaaaaaaaULL
        }};
        const cc512i a = {{
            0, 1, mask52, mask52 - 1, u64{1} << 51,
            std::numeric_limits<u64>::max(), (u64{0xfff} << 52) | 1, (u64{1} << 52) | mask52
        }};
        const cc512i b = {{
            mask52, mask52, mask52, u64{1} << 51, u64{1} << 26,
            std::numeric_limits<u64>::max(), (u64{0xabc} << 52) | mask52, (u64{1} << 63) | 1
        }};
        check_madd52_case(acc, a, b);
    }

    TEST_F(Avx512CompatTest, Madd52Random)
    {
        for (int iteration = 0; iteration < 256; ++iteration) {
            check_madd52_case(random_cc512(), random_cc512(), random_cc512());
        }
    }
} // namespace blst::test
