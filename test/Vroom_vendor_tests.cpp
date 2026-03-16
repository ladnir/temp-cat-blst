#include <gtest/gtest.h>

#include <array>
#include <sstream>

#include "InsecureRng.h"
#include "blst/crns_batch/VroomLiteral.h"

#if defined(BLST_CRNS_BATCH_VROOM_VENDOR)

using namespace blst;

namespace
{
std::string fp_hex(const Fp& value)
{
    Fp canonical{};
    Fp::from_mont(canonical, value);
    std::ostringstream oss;
    oss << std::hex;
    for (int i = static_cast<int>(Fp::num_limb()) - 1; i >= 0; --i) {
        oss << canonical[static_cast<size_t>(i)] << ":";
    }
    return oss.str();
}

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
} // namespace

TEST(VroomVendor, BatchModMulMatchesScalar)
{
    InsecureRng rng(0x10293847u);
    std::array<Fp, 8> a{};
    std::array<Fp, 8> b{};
    std::array<Fp, 8> got{};
    std::array<Fp, 8> expected{};

    for (size_t lane = 0; lane < 8; ++lane) {
        rand_redc_fp(a[lane], rng);
        rand_redc_fp(b[lane], rng);
        Fp::mul(expected[lane], a[lane], b[lane]);
    }

    blst::crns_batch::vroom_literal::VroomFpPack8 pa{};
    blst::crns_batch::vroom_literal::VroomFpPack8 pb{};
    blst::crns_batch::vroom_literal::VroomFpPack8 pr{};
    blst::crns_batch::vroom_literal::from_aos(pa, a.data());
    blst::crns_batch::vroom_literal::from_aos(pb, b.data());
    blst::crns_batch::vroom_literal::mul(pr, pa, pb);
    blst::crns_batch::vroom_literal::to_aos(got.data(), pr);

    for (size_t lane = 0; lane < 8; ++lane) {
        EXPECT_TRUE(Fp::equal(got[lane], expected[lane])) << "lane=" << lane;
    }
}

TEST(VroomVendor, Fp12MulMatchesScalar)
{
    InsecureRng rng(0x56473829u);
    for (int round = 0; round < 8; ++round) {
        Fp12 a{};
        Fp12 b{};
        Fp12 expected{};
        rand_redc_fp12(a, rng);
        rand_redc_fp12(b, rng);

        Fp12::mul(expected, a, b);

        const auto va = blst::crns_batch::vroom_literal::fp12_from_blst(a);
        const auto vb = blst::crns_batch::vroom_literal::fp12_from_blst(b);
        const auto vr = blst::crns_batch::vroom_literal::fp12_mul(va, vb);
        const auto got = blst::crns_batch::vroom_literal::fp12_to_blst(vr);

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 2; ++k) {
                    ASSERT_TRUE(Fp::equal(got[i][j][k], expected[i][j][k]))
                        << "round=" << round
                        << " coeff=" << i << "/" << j << "/" << k
                        << " got=" << fp_hex(got[i][j][k])
                        << " expected=" << fp_hex(expected[i][j][k]);
                }
            }
        }
    }
}

#endif
