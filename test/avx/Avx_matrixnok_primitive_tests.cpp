#include "blst/avx/generated_matrixnok_constants.hpp"
#include "blst/avx/cpu/vector/conversion.hpp"
#include "blst/avx/rns/elementwise_reduce.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>

namespace blst::avx::test
{
    namespace matrixnok_constants = matrixnok::generated;

    static_assert(ConvertToRNS<381, 8>::indigits == 8);
    static_assert(ConvertToRNS<381, 8, 50>::indigits == 8);
    static_assert(ConvertToRNS<401, 8>::indigits == 8);
    static_assert(ConvertToRNS<401, 8, 50>::indigits == 9);
    static_assert(matrixnok_constants::BLS12_381_RNS_DIM == 8);
    static_assert(matrixnok_constants::BLS12_381_CONVERT_TO_RADIX_BITS == 50);
    static_assert(matrixnok_constants::moduli1.size() == 8);
    static_assert(matrixnok_constants::moduli2.size() == 8);
    static_assert(matrixnok_constants::r1_perm_mat_fixed.size() == 8);

    TEST(AvxMatrixNoKPrimitives, LowLaneExtractionAndBroadcast)
    {
        const std::array<uint64_t, 8> low = {
            (uint64_t{7} << 52) | 5, 11, 13, 17, 19, 23, 29, 31};
        const std::array<uint64_t, 8> high = {
            3, 37, 41, 43, 47, 53, 59, 61};

        const ShortVector low_short = ShortVector::from_avx_lo(AVXVector<8>(low));
        const ShortVector high_short = ShortVector::from_avx_lo(AVXVector<8>(high));
        const auto broadcast = high_short.add(low_short.srli(52)).broadcast_to_avx<8>().to_array();

        for (const uint64_t lane : broadcast) {
            EXPECT_EQ(lane, 10);
        }
        EXPECT_EQ(AVXVector<8>(low).extract0(), low[0]);
    }

    template<int64_t upper_bound, int64_t stop, int log_multiples>
    void check_min_sub_schedule()
    {
        constexpr int limbs = 8;
        constexpr int element_bits = 50;
        const std::array<uint64_t, limbs> moduli = {
            (uint64_t{1} << element_bits) - 1,
            (uint64_t{1} << element_bits) - 3,
            (uint64_t{1} << element_bits) - 5,
            (uint64_t{1} << element_bits) - 7,
            (uint64_t{1} << element_bits) - 11,
            (uint64_t{1} << element_bits) - 13,
            (uint64_t{1} << element_bits) - 17,
            (uint64_t{1} << element_bits) - 19,
        };
        const MontgomeryReduce2<limbs, element_bits, log_multiples> reducer(moduli);

        constexpr std::array<int64_t, 16> quotient_seeds = {
            0, 1, 2, 3, 4, 5, 7, 8, 15, 16, 31, 32,
            upper_bound / 2, upper_bound - 3, upper_bound - 2, upper_bound - 1,
        };
        for (const int64_t seed : quotient_seeds) {
            const int64_t quotient = seed < 0 ? 0 : (seed >= upper_bound ? upper_bound - 1 : seed);
            std::array<uint64_t, limbs> values{};
            for (int lane = 0; lane < limbs; ++lane) {
                const uint64_t remainder = (lane & 1) == 0 ? 0 : moduli[lane] - 1;
                values[lane] = static_cast<uint64_t>(quotient) * moduli[lane] + remainder;
            }

            using Input = BoundedElement<Bounds<0, upper_bound>, limbs, element_bits>;
            const auto reduced = reducer.template min_sub_reduce<upper_bound, stop>(Input(AVXVector<limbs>(values)));
            const auto actual = reduced.data.to_array();
            for (int lane = 0; lane < limbs; ++lane) {
                EXPECT_LT(actual[lane], static_cast<uint64_t>(stop) * moduli[lane]);
                EXPECT_EQ(actual[lane] % moduli[lane], values[lane] % moduli[lane]);
            }
        }
    }

    TEST(AvxMatrixNoKPrimitives, ConditionalSubtractionCoversConfiguredGap)
    {
        check_min_sub_schedule<2, 1, 1>();
        check_min_sub_schedule<4, 1, 2>();
        check_min_sub_schedule<40, 1, 6>();
        check_min_sub_schedule<100, 95, 3>();
    }
} // namespace blst::avx::test
