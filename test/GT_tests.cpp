#include "gtest/gtest.h"

#include "blst/curve/e1.h"
#include "blst/curve/e2.h"
#include "blst/curve/pairing.h"
#include "blst/fields/Fp12.h"
#include "blst/sign/aggregate.h"

using namespace blst;

class GTTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize points
        g1_point = G1Point::generator();
        g2_point = G2Point::generator();

        // Convert to affine
        G1Point::to_affine(g1_affine, g1_point);
        G2Point::to_affine(g2_affine, g2_point);

        // Calculate pairing
        blst_miller_loop(gt_element, &g2_affine, &g1_affine);
        Fp12::final_exp(gt_element, gt_element);
    }

    G1Point g1_point;
    G2Point g2_point;
    G1AffinePoint g1_affine;
    G2AffinePoint g2_affine;
    Fp12 gt_element;  // An element in GT
    Fp12 temp;        // For temporary results
};

TEST_F(GTTest, Identity) {
    // Test identity element in GT
    auto identity = Fp12::one();

    EXPECT_TRUE(Fp12::is_one(identity))
        << "Identity element should be 1";

    // Multiplying any element by identity should leave it unchanged
    Fp12::mul(temp, gt_element, identity);
    EXPECT_TRUE(Fp12::equal(temp, gt_element))
        << "a * 1 should equal a in GT";
}

TEST_F(GTTest, Multiplication) {
    // Create another element for testing
    Fp12 another_gt;
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    
    blst_miller_loop(another_gt, &g2_affine, &g1_double_affine);
    Fp12::final_exp(another_gt, another_gt);

    // Test multiplication
    Fp12 product;
    Fp12::mul(product, gt_element, another_gt);

    // The product should not be equal to either of the original elements
    EXPECT_FALSE(Fp12::equal(product, gt_element))
        << "Product should not equal first element";
    EXPECT_FALSE(Fp12::equal(product, another_gt))
        << "Product should not equal second element";

    // Test commutativity: a * b = b * a
    Fp12 product2;
    Fp12::mul(product2, another_gt, gt_element);
    EXPECT_TRUE(Fp12::equal(product, product2))
        << "Multiplication should be commutative in GT";

    // Test associativity: (a * b) * c = a * (b * c)
    // Create a third element
    Fp12 third_gt;
    G1Point g1_triple;
    G1Point::add(g1_triple, g1_double, g1_point);

    G1AffinePoint g1_triple_affine;
    G1Point::to_affine(g1_triple_affine, g1_triple);
    
    blst_miller_loop(third_gt, &g2_affine, &g1_triple_affine);
    Fp12::final_exp(third_gt, third_gt);

    Fp12 temp1, temp2;
    // (a * b) * c
    Fp12::mul(temp1, gt_element, another_gt);
    Fp12::mul(temp1, temp1, third_gt);

    // a * (b * c)
    Fp12::mul(temp2, another_gt, third_gt);
    Fp12::mul(temp2, gt_element, temp2);

    EXPECT_TRUE(Fp12::equal(temp1, temp2))
        << "Multiplication should be associative in GT";
}

TEST_F(GTTest, Inversion) {
    // Test inversion (a * a⁻¹ = 1)
    Fp12 inverse;
    Fp12::inverse(inverse, gt_element);

    Fp12 product;
    Fp12::mul(product, gt_element, inverse);

    auto identity = Fp12::one();

    EXPECT_TRUE(Fp12::equal(product, identity))
        << "Element multiplied by its inverse should equal identity";
}

TEST_F(GTTest, Squaring) {
    // Test that squaring equals multiplying by self
    Fp12 square;
    Fp12::sqr(square, gt_element);

    Fp12 product;
    Fp12::mul(product, gt_element, gt_element);

    EXPECT_TRUE(Fp12::equal(square, product))
        << "a² should equal a * a";

    // Test multiple squaring (to test repeated squaring optimization if available)
    Fp12 repeated_square;
    Fp12::sqr(repeated_square, gt_element);
    Fp12::sqr(repeated_square, repeated_square);

    // This should equal gt_element⁴
    Fp12 pow4;
    Fp12::mul(pow4, gt_element, gt_element);
    Fp12::mul(pow4, pow4, pow4);

    EXPECT_TRUE(Fp12::equal(repeated_square, pow4))
        << "((a²)²) should equal a⁴";
}

TEST_F(GTTest, Exponentiation) {
    // Test exponentiation using small exponents
    Fp12 result, expected;

    // Test e^0 = 1
    limb_t zero[4] = { 0, 0, 0, 0 };
    Fp12::exp(result, gt_element, zero, sizeof(zero) * 8);
    Fp12 expected_identity = { 0 };
	expected_identity = Fp12::one();
    EXPECT_TRUE(Fp12::equal(result, expected_identity))
        << "e^0 should equal 1";

    // Test e^1 = e
    limb_t one[4] = { 1, 0, 0, 0 };
    Fp12::exp(result, gt_element, one, sizeof(one) * 8);
    EXPECT_TRUE(Fp12::equal(result, gt_element))
        << "e^1 should equal e";

    // Test e^2 = e * e
    limb_t two[4] = { 2, 0, 0, 0 };
    Fp12::exp(result, gt_element, two, sizeof(two) * 8);
    Fp12::sqr(expected, gt_element);
    EXPECT_TRUE(Fp12::equal(result, expected))
        << "e^2 should equal e * e";

    // Test e^3 = e * e * e
    limb_t three[4] = { 3, 0, 0, 0 };
    Fp12::exp(result, gt_element, three, sizeof(three) * 8);
    Fp12::mul(expected, expected, gt_element);
    EXPECT_TRUE(Fp12::equal(result, expected))
        << "e^3 should equal e * e * e";
}

//TEST_F(GTTest, Serialization) {
//    // Test serialization/deserialization of GT elements
//    Fp12 original = { 0 };  // Copy the gt_element
//    memcpy(original, gt_element, sizeof(Fp12));
//
//    // Serialize to bytes
//    unsigned char serialized[12 * 48];  // Fp12 consists of 12 base field elements
//    blst_fp12_serialize(serialized, original);
//
//    // Deserialize back
//    Fp12 deserialized;
//    EXPECT_TRUE(blst_fp12_deserialize(deserialized, serialized))
//        << "Deserialization should succeed";
//
//    // The original and deserialized values should be equal
//    EXPECT_TRUE(Fp12::equal(original, deserialized))
//        << "Deserialized value should equal original";
//}

TEST_F(GTTest, BilinearityProperty) {
    // Test the bilinearity property: e(a*P, Q) = e(P, a*Q) = e(P,Q)^a

    // Create scaled points
    G1Point g1_double;
    G2Point g2_double;
    G1Point::doubling(g1_double, g1_point);
    G2Point::doubling(g2_double, g2_point);

    G1AffinePoint g1_double_affine;
    G2AffinePoint g2_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_double_affine, g2_double);

    // e(2P, Q)
    Fp12 e_2P_Q;
    blst_miller_loop(e_2P_Q, &g2_affine, &g1_double_affine);
    Fp12::final_exp(e_2P_Q, e_2P_Q);

    // e(P, 2Q)
    Fp12 e_P_2Q;
    blst_miller_loop(e_P_2Q, &g2_double_affine, &g1_affine);
    Fp12::final_exp(e_P_2Q, e_P_2Q);

    // e(P,Q)^2
    Fp12 e_PQ_squared;
    Fp12::sqr(e_PQ_squared, gt_element);

    // All three should be equal due to bilinearity
    EXPECT_TRUE(Fp12::equal(e_2P_Q, e_P_2Q))
        << "e(2P, Q) should equal e(P, 2Q)";
    EXPECT_TRUE(Fp12::equal(e_2P_Q, e_PQ_squared))
        << "e(2P, Q) should equal e(P,Q)^2";
}



// NEW TESTS BELOW

TEST_F(GTTest, GroupMembership) {
    // Test that the pairing result is in the correct group
    EXPECT_TRUE(Fp12::in_group(gt_element))
        << "Pairing result should be in the GT group";

    // Identity should be in the group
    Fp12 identity;
    identity = Fp12::one();
    EXPECT_TRUE(Fp12::in_group(identity))
        << "Identity element should be in the GT group";

    // Test that a non-GT element (before final_exp) is not in the group
    Fp12 not_finalized;
    blst_miller_loop(not_finalized, &g2_affine, &g1_affine);
    EXPECT_FALSE(Fp12::in_group(not_finalized))
        << "Miller loop result without final_exp should not be in GT group";

    // Squaring, multiplication, inversion should preserve group membership
    Fp12 squared;
    Fp12::sqr(squared, gt_element);
    EXPECT_TRUE(Fp12::in_group(squared))
        << "Squaring should preserve group membership";

    Fp12 inverse;
    Fp12::inverse(inverse, gt_element);
    EXPECT_TRUE(Fp12::in_group(inverse))
        << "Inversion should preserve group membership";
}

TEST_F(GTTest, MultiPairing) {
    // Set up multiple points for testing
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_points[2];
    G1Point::to_affine(g1_points[0], g1_point);
    G1Point::to_affine(g1_points[1], g1_double);

    G2AffinePoint g2_points[2];
    G2Point::to_affine(g2_points[0], g2_point);
    G2Point::to_affine(g2_points[1], g2_point);

    // Create arrays of pointers as required by blst_miller_loop_n
    const G1AffinePoint* g1_ptrs[3] = { &g1_points[0], &g1_points[1], nullptr };
    const G2AffinePoint* g2_ptrs[3] = { &g2_points[0], &g2_points[1], nullptr };

    // Calculate e(P1,Q1) and e(P2,Q2) individually
    Fp12 e1, e2, product, multi;
    blst_miller_loop(e1, &g2_points[0], &g1_points[0]);
    blst_miller_loop(e2, &g2_points[1], &g1_points[1]);

    // Calculate their product
    Fp12::mul(product, e1, e2);
    Fp12::final_exp(product, product);

    // Calculate using multi-pairing
    blst_miller_loop_n(multi, g2_ptrs, g1_ptrs, 2);
    Fp12::final_exp(multi, multi);

    // They should match
    EXPECT_TRUE(Fp12::equal(product, multi))
        << "Multi-pairing result should match individual pairings";
}

TEST_F(GTTest, FinalVerify) {
    // Test finalverify function (equality after final exponentiation)

    // Create two different pairings that should be equal after final_exp
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    // Create a "half" point in G2 by multiplying by inverse of 2 (if available)
    // For test purposes, we'll just use a different point
    G2Point g2_another;
    G2Point::add(g2_another, g2_point, g2_point);  // Just a different point

    G1AffinePoint g1_double_affine;
    G2AffinePoint g2_half_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_half_affine, g2_another);

    // Calculate two different pairings
    Fp12 pair1, pair2;
    blst_miller_loop(pair1, &g2_affine, &g1_double_affine);
    blst_miller_loop(pair2, &g2_half_affine, &g1_affine);

    // They should not be equal before final_exp
    EXPECT_FALSE(Fp12::equal(pair1, pair2))
        << "Different pairings should not be equal before final_exp";

    // Apply final exponentiation
    Fp12 final1, final2;
    Fp12::final_exp(final1, pair1);
    Fp12::final_exp(final2, pair2);

    // Direct equality check
    bool direct_equal = Fp12::equal(final1, final2);

    // Using finalverify
    bool finalverify_equal = blst_fp12_finalverify(final1, final2);

    // Verify that finalverify gives same result as equality check
    // (in general they might differ due to optimizations, but for our simple test they should match)
    EXPECT_EQ(direct_equal, finalverify_equal)
        << "finalverify should agree with is_equal in basic cases";
}

TEST_F(GTTest, PairingWithPrecomputation) {
    // Test using precomputed lines for pairing
    Fp6 lines[68];  // Size from the API

    // Precompute lines for pairing with g2_point
    blst_precompute_lines(lines, g2_affine);

    // Compute pairing using precomputed lines
    Fp12 pre_result;
    blst_miller_loop_lines(pre_result, lines, g1_affine);
    Fp12::final_exp(pre_result, pre_result);

    // Verify against standard pairing computation
    EXPECT_TRUE(Fp12::equal(gt_element, pre_result))
        << "Pairing with precomputation should match direct computation";

    // Verify that we can reuse the precomputed lines with a different point
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    
    Fp12 direct_double_result;
    blst_miller_loop(direct_double_result, &g2_affine, &g1_double_affine);
    Fp12::final_exp(direct_double_result, direct_double_result);

    Fp12 pre_double_result;
    blst_miller_loop_lines(pre_double_result, lines, g1_double_affine);
    Fp12::final_exp(pre_double_result, pre_double_result);

    EXPECT_TRUE(Fp12::equal(direct_double_result, pre_double_result))
        << "Precomputed pairing should work with different G1 points";
}

TEST_F(GTTest, MillerLoopOnly) {
    // Test Miller loop without final exponentiation
    Fp12 miller_result1, miller_result2;
    
    // Compute Miller loop for same pairing twice
    blst_miller_loop(miller_result1, &g2_affine, &g1_affine);
    blst_miller_loop(miller_result2, &g2_affine, &g1_affine);
    
    // Should be equal
    EXPECT_TRUE(Fp12::equal(miller_result1, miller_result2))
        << "Miller loop should be deterministic";
    
    // Should not be in GT group before final exp
    EXPECT_FALSE(Fp12::in_group(miller_result1))
        << "Miller loop result should not be in GT before final exponentiation";
    
    // Apply final exp and verify it's now in GT
    Fp12 final_result;
    Fp12::final_exp(final_result, miller_result1);
    EXPECT_TRUE(Fp12::in_group(final_result))
        << "Result should be in GT after final exponentiation";
}

TEST_F(GTTest, PairingIdentityProperties) {
    // Test pairing with identity elements
    
    // Create identity point in G1 (point at infinity)
    G1Point g1_identity;
    memset(&g1_identity, 0, sizeof(g1_identity));
    // Set z coordinate to 0 to make it point at infinity
    //g1_identity.z[0] = 0;
    
    G1AffinePoint g1_identity_affine;
    G1Point::to_affine(g1_identity_affine, g1_identity);
    
    // Pairing with identity should give identity in GT
    Fp12 pairing_with_identity;
    blst_miller_loop(pairing_with_identity, &g2_affine, &g1_identity_affine);
    Fp12::final_exp(pairing_with_identity, pairing_with_identity);
    
    EXPECT_TRUE(Fp12::is_one(pairing_with_identity))
        << "Pairing with identity point should give identity in GT";
}