#define BLST_AVX_SKIP_TEST_VALUES

#include <gtest/gtest.h>

#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/point.h"
#include "avx_test_utils.h"

#include <array>
#include <cstring>

using namespace blst_avx_test;

namespace {

TEST(AvxWrapper, G1GeneratorSerializationMatchesScalar)
{
    unsigned char scalar_bytes[96] = {};
    unsigned char avx_bytes[96] = {};

    blst::G1AffinePoint::serialize(scalar_bytes, blst::G1AffinePoint::generator());
    blst::avx::G1AffinePoint::serialize(avx_bytes, blst::avx::G1AffinePoint::generator());

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G1DoublingMatchesScalar)
{
    blst::G1Point scalar_dbl = {};
    blst::G1Point::doubling(scalar_dbl, blst::G1Point::generator());

    blst::avx::G1Point avx_dbl = {};
    blst::avx::G1Point::doubling(avx_dbl, blst::avx::G1Point::generator());

    unsigned char scalar_bytes[96] = {};
    unsigned char avx_bytes[96] = {};
    blst::G1Point::serialize(scalar_bytes, scalar_dbl);
    blst::avx::G1Point::serialize(avx_bytes, avx_dbl);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G1ScalarMultMatchesScalar)
{
    const std::array<blst::byte, 32> scalar = {1, 2, 3, 4, 5, 6, 7, 8};

    blst::G1Point scalar_out = {};
    blst::G1Point::mult(scalar_out, blst::G1Point::generator(), scalar.data(), 256);

    blst::avx::G1Point avx_out = {};
    blst::avx::G1Point::mult(avx_out, blst::avx::G1Point::generator(), scalar.data(), 256);

    unsigned char scalar_bytes[96] = {};
    unsigned char avx_bytes[96] = {};
    blst::G1Point::serialize(scalar_bytes, scalar_out);
    blst::avx::G1Point::serialize(avx_bytes, avx_out);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G2GeneratorSerializationMatchesScalar)
{
    unsigned char scalar_bytes[192] = {};
    unsigned char avx_bytes[192] = {};

    blst::G2AffinePoint::serialize(scalar_bytes, blst::G2AffinePoint::generator());
    blst::avx::G2AffinePoint::serialize(avx_bytes, blst::avx::G2AffinePoint::generator());

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G2DoublingMatchesScalar)
{
    blst::G2Point scalar_dbl = {};
    blst::G2Point::doubling(scalar_dbl, blst::G2Point::generator());

    blst::avx::G2Point avx_dbl = {};
    blst::avx::G2Point::doubling(avx_dbl, blst::avx::G2Point::generator());

    unsigned char scalar_bytes[192] = {};
    unsigned char avx_bytes[192] = {};
    blst::G2Point::serialize(scalar_bytes, scalar_dbl);
    blst::avx::G2Point::serialize(avx_bytes, avx_dbl);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G2ScalarMultMatchesScalar)
{
    const std::array<blst::byte, 32> scalar = {9, 8, 7, 6, 5, 4, 3, 2};

    blst::G2Point scalar_out = {};
    blst::G2Point::mult(scalar_out, blst::G2Point::generator(), scalar.data(), 256);

    blst::avx::G2Point avx_out = {};
    blst::avx::G2Point::mult(avx_out, blst::avx::G2Point::generator(), scalar.data(), 256);

    unsigned char scalar_bytes[192] = {};
    unsigned char avx_bytes[192] = {};
    blst::G2Point::serialize(scalar_bytes, scalar_out);
    blst::avx::G2Point::serialize(avx_bytes, avx_out);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTRoundTripSerializationMatchesScalar)
{
    const blst::GT scalar_gt = make_scalar_final_pairing(
        make_scalar_g1_generator_affine(),
        make_scalar_g2_generator_affine());

    unsigned char scalar_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_gt);

    blst::avx::GT avx_gt = {};
    blst::avx::GT::deserialize(avx_gt, scalar_bytes);

    unsigned char avx_bytes[12 * 48] = {};
    blst::avx::GT::serialize(avx_bytes, avx_gt);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTOperationsMatchScalar)
{
    const blst::GT scalar_a = make_scalar_final_pairing(
        make_scalar_g1_generator_affine(),
        make_scalar_g2_generator_affine());
    const blst::GT scalar_b = make_scalar_final_pairing(
        make_scalar_g1_doubled_affine(),
        make_scalar_g2_generator_affine());

    unsigned char scalar_a_bytes[12 * 48] = {};
    unsigned char scalar_b_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_a_bytes, scalar_a);
    blst::GT::serialize(scalar_b_bytes, scalar_b);

    blst::avx::GT avx_a = {};
    blst::avx::GT avx_b = {};
    blst::avx::GT::deserialize(avx_a, scalar_a_bytes);
    blst::avx::GT::deserialize(avx_b, scalar_b_bytes);

    blst::GT scalar_mul = {};
    blst::GT::add(scalar_mul, scalar_a, scalar_b);

    blst::avx::GT avx_mul = {};
    blst::avx::GT::add(avx_mul, avx_a, avx_b);

    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_mul);
    blst::avx::GT::serialize(avx_bytes, avx_mul);
    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));

    blst::GT scalar_sq = {};
    blst::GT::doubling(scalar_sq, scalar_a);

    blst::avx::GT avx_sq = {};
    blst::avx::GT::doubling(avx_sq, avx_a);

    blst::GT::serialize(scalar_bytes, scalar_sq);
    blst::avx::GT::serialize(avx_bytes, avx_sq);
    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTScalarMultMatchesScalar)
{
    const blst::GT scalar_a = make_scalar_final_pairing(
        make_scalar_g1_generator_affine(),
        make_scalar_g2_generator_affine());
    const std::array<blst::byte, 32> scalar = {1, 0, 2, 0, 3, 0, 4, 0};

    unsigned char scalar_a_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_a_bytes, scalar_a);

    blst::avx::GT avx_a = {};
    blst::avx::GT::deserialize(avx_a, scalar_a_bytes);

    blst::GT scalar_out = {};
    blst::GT::mult(scalar_out, scalar_a, scalar.data(), 256);

    blst::avx::GT avx_out = {};
    blst::avx::GT::mult(avx_out, avx_a, scalar.data(), 256);

    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_out);
    blst::avx::GT::serialize(avx_bytes, avx_out);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

}  // namespace
