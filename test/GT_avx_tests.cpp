#include <gtest/gtest.h>

#include <array>

#include "InsecureRng.h"
#include "blst/curve/GT_avx.h"

using namespace blst;

namespace
{
void rand_redc_fp(Fp& v, InsecureRng& rng)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = rng.next_u64();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

void rand_redc_fp2(Fp2& v, InsecureRng& rng)
{
    rand_redc_fp(v[0], rng);
    rand_redc_fp(v[1], rng);
}

void rand_redc_fp12(Fp12& v, InsecureRng& rng)
{
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rand_redc_fp2(v[i][j], rng);
        }
    }
}

void make_random_gt(GT& out, InsecureRng& rng)
{
    Fp12 seed{};
    rand_redc_fp12(seed, rng);
    Fp12::final_exp(out.data, seed);
}
} // namespace

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

TEST(GTBatchAvx, PackRoundTripAndSelect)
{
    InsecureRng rng;
    std::array<GT, 8> a{};
    std::array<GT, 8> b{};
    std::array<GT, 8> out{};
    GTPack8x1 pa{};
    GTPack8x1 pb{};
    GTPack8x1 pr{};

    for (size_t lane = 0; lane < 8; ++lane) {
        make_random_gt(a[lane], rng);
        make_random_gt(b[lane], rng);
    }

    GTPack8x1::from_aos(pa, a.data());
    GTPack8x1::to_aos(out.data(), pa);
    for (size_t lane = 0; lane < 8; ++lane) {
        EXPECT_TRUE(GT::equal(out[lane], a[lane]));
    }

    GTPack8x1::from_aos(pb, b.data());
    const std::array<u8, 5> masks = {0x00, 0x55, 0xAA, 0xFF, 0x93};
    for (const u8 mask : masks) {
        GTPack8x1::select(pr, pa, pb, mask);
        GTPack8x1::to_aos(out.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            const GT& expected = ((mask >> lane) & 1u) ? b[lane] : a[lane];
            EXPECT_TRUE(GT::equal(out[lane], expected));
        }
    }
}

TEST(GTBatchAvx, PackPartialRoundTrip)
{
    InsecureRng rng;
    const GT fill = GT::identity();

    for (size_t n_active = 0; n_active <= 8; ++n_active) {
        std::array<GT, 8> in{};
        std::array<GT, 8> out{};
        std::array<GT, 8> full{};
        GTPack8x1 pack{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(in[lane], rng);
        }

        GTPack8x1::from_aos_partial(pack, in.data(), n_active, fill);
        GTPack8x1::to_aos_partial(out.data(), n_active, pack);
        GTPack8x1::to_aos(full.data(), pack);

        for (size_t lane = 0; lane < n_active; ++lane) {
            EXPECT_TRUE(GT::equal(out[lane], in[lane]));
            EXPECT_TRUE(GT::equal(full[lane], in[lane]));
        }
        for (size_t lane = n_active; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(full[lane], fill));
        }
    }
}

TEST(GTBatchAvx, AddMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 48; ++iter) {
        std::array<GT, 8> a{};
        std::array<GT, 8> b{};
        std::array<GT, 8> got{};
        std::array<GT, 8> expected{};
        GTPack8x1 pa{};
        GTPack8x1 pb{};
        GTPack8x1 pr{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(a[lane], rng);
            make_random_gt(b[lane], rng);
            GT::add(expected[lane], a[lane], b[lane]);
        }

        GTPack8x1::from_aos(pa, a.data());
        GTPack8x1::from_aos(pb, b.data());

        GTPack8x1::add(pr, pa, pb);
        GTPack8x1::to_aos(got.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }

        GTPack8x1::add(pr, pa, pb);
        GTPack8x1::to_aos(got.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

TEST(GTBatchAvx, AddAliasing)
{
    InsecureRng rng;

    for (int iter = 0; iter < 48; ++iter) {
        std::array<GT, 8> a{};
        std::array<GT, 8> b{};
        std::array<GT, 8> expected{};
        std::array<GT, 8> got{};
        GTPack8x1 pa{};
        GTPack8x1 pb{};
        GTPack8x1 pa_orig{};
        GTPack8x1 pb_orig{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(a[lane], rng);
            make_random_gt(b[lane], rng);
            GT::add(expected[lane], a[lane], b[lane]);
        }

        GTPack8x1::from_aos(pa, a.data());
        GTPack8x1::from_aos(pb, b.data());
        pa_orig = pa;
        pb_orig = pb;

        GTPack8x1::add(pa, pa, pb);
        GTPack8x1::to_aos(got.data(), pa);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }

        GTPack8x1::add(pb, pa_orig, pb);
        GTPack8x1::to_aos(got.data(), pb);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

TEST(GTBatchAvx, DoublingMatchesScalarAndAliasing)
{
    InsecureRng rng;

    for (int iter = 0; iter < 48; ++iter) {
        std::array<GT, 8> a{};
        std::array<GT, 8> got{};
        std::array<GT, 8> expected{};
        GTPack8x1 pa{};
        GTPack8x1 pinout{};
        GTPack8x1 pr{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(a[lane], rng);
            ASSERT_TRUE(GT::is_cyclotomic(a[lane]));
            GT::doubling(expected[lane], a[lane]);
        }

        GTPack8x1::from_aos(pa, a.data());

        GTPack8x1::doubling(pr, pa);
        GTPack8x1::to_aos(got.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }

        pinout = pa;
        GTPack8x1::doubling(pinout, pinout);
        GTPack8x1::to_aos(got.data(), pinout);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

#else

TEST(GTBatchAvx, Disabled)
{
    SUCCEED();
}

#endif
