#include "gtest/gtest.h"

#include "blst/curve/point.h"
#include "blst/curve/e1.h"
#include "blst/curve/e2.h"
#include "blst/curve/pairing.h"
using namespace blst;


TEST(PairingTest, BasicPairing) {
    // Test basic pairing functionality
    G1Point g1_point = G1Point::generator();
    G2Point g2_point = G2Point::generator();

    // Create affine versions of the points
    G1AffinePoint g1_affine;
    G2AffinePoint g2_affine;
    G1Point::to_affine(g1_affine, g1_point);
    G2Point::to_affine(g2_affine, g2_point);

    // Compute pairing e(G1, G2)
    Fp12 pairing_result;
    blst_miller_loop(pairing_result, &g2_affine, &g1_affine);
    Fp12::final_exp(pairing_result, pairing_result);

    // The pairing should not be equal to 1 (non-degeneracy)
    EXPECT_FALSE(Fp12::is_one(pairing_result))
        << "Pairing should be non-degenerate";

    // Test bilinearity: e(aP, Q) = e(P, aQ) = e(P,Q)^a
    // Let's use a = 2 for simplicity
    G1Point g1_double;
    G2Point g2_double;
    G1Point::doubling(g1_double, g1_point);
    G2Point::doubling(g2_double, g2_point);

    // Convert to affine
    G1AffinePoint g1_double_affine;
    G2AffinePoint g2_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_double_affine, g2_double);

    // Compute e(2*G1, G2)
    Fp12 pairing_2g1_g2;
    blst_miller_loop(pairing_2g1_g2, &g2_affine, &g1_double_affine);
    Fp12::final_exp(pairing_2g1_g2, pairing_2g1_g2);

    // Compute e(G1, 2*G2)
    Fp12 pairing_g1_2g2;
    blst_miller_loop(pairing_g1_2g2, &g2_double_affine, &g1_affine);
    Fp12::final_exp(pairing_g1_2g2, pairing_g1_2g2);

    // Compute e(G1, G2)^2
    Fp12 pairing_squared;
    Fp12::sqr(pairing_squared, pairing_result);

    // All three results should be equal
    EXPECT_TRUE(Fp12::equal(pairing_2g1_g2, pairing_g1_2g2))
        << "e(2*G1, G2) != e(G1, 2*G2)";

    EXPECT_TRUE(Fp12::equal(pairing_2g1_g2, pairing_squared))
        << "e(2*G1, G2) != e(G1, G2)^2";

    // Test pairing with identity points
    G1Point identity_e1;
    G2Point identity_e2;
    memset(&identity_e1, 0, sizeof(identity_e1));
    memset(&identity_e2, 0, sizeof(identity_e2));

    G1AffinePoint identity_e1_affine;
    G2AffinePoint identity_e2_affine;
    G1Point::to_affine(identity_e1_affine, identity_e1);
    G2Point::to_affine(identity_e2_affine, identity_e2);

    // e(G1, O) should equal 1
    Fp12 pairing_with_identity;
    blst_miller_loop(pairing_with_identity, &identity_e2_affine, &g1_affine);
    Fp12::final_exp(pairing_with_identity, pairing_with_identity);

    EXPECT_TRUE(Fp12::is_one(pairing_with_identity))
        << "e(G1, O) should equal 1";

    // e(O, G2) should equal 1
    blst_miller_loop(pairing_with_identity,& g2_affine, &identity_e1_affine);
    Fp12::final_exp(pairing_with_identity, pairing_with_identity);

    EXPECT_TRUE(Fp12::is_one(pairing_with_identity))
        << "e(O, G2) should equal 1";
}

TEST(PairingTest, MultiPairing) {
    // Test multi-pairing functionality

    // Create some sample points
    G1Point g1 = G1Point::generator();
    G2Point g2 = G2Point::generator();

    G1Point g1_double;
    G2Point g2_double;
    G1Point::doubling(g1_double, g1);
    G2Point::doubling(g2_double, g2);

    // Convert to affine
    G1AffinePoint g1_affine, g1_double_affine;
    G2AffinePoint g2_affine, g2_double_affine;
    G1Point::to_affine(g1_affine, g1);
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_affine, g2);
    G2Point::to_affine(g2_double_affine, g2_double);

    // Test that e(P1,Q1)·e(P2,Q2) = e(P1+P2,Q1+Q2)·e(P1,Q2)⁻¹·e(P2,Q1)⁻¹

    // Compute e(G1, G2)
    Fp12 e_g1_g2;
    blst_miller_loop(e_g1_g2, &g2_affine, &g1_affine);
    Fp12::final_exp(e_g1_g2, e_g1_g2);

    // Compute e(2*G1, 2*G2)
    Fp12 e_2g1_2g2;
    blst_miller_loop(e_2g1_2g2, &g2_double_affine,&g1_double_affine);
    Fp12::final_exp(e_2g1_2g2, e_2g1_2g2);

    // Compute e(G1, G2) * e(2*G1, 2*G2) using Fp12_mul
    Fp12 product;
    Fp12::mul(product, e_g1_g2, e_2g1_2g2);

    // Test multi-pairing API if available
    G1AffinePoint* g1_points[2] = { &g1_affine, &g1_double_affine };
    G2AffinePoint* g2_points[2] = { &g2_affine, &g2_double_affine };

    Fp12 multi_result;
    blst_miller_loop_n(multi_result, g2_points, g1_points, 2);
    Fp12::final_exp(multi_result, multi_result);

    // The results should be equal
    EXPECT_TRUE(Fp12::equal(product, multi_result))
        << "Multi-pairing result does not match individual pairings";
}