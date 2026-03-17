#include "../cpu/precompute/gmp_wrapper.hpp"
extern "C" {
#include "../blst/vect.h"
#include "../blst/fields.h"
#include "../blst/consts.h"
#include "../blst/point.h"
#include "../blst/bytes.h"
#include "../blst/ec_mult.h"
}
#include "scalar_mult.hpp"
#include "bounded_ring.hpp"
#include "fp2.hpp"
#include "conversion_inversion.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <cstring>
#include <vector>

// BLS12-381 modulus q
const char* bls12_381_modulus_hex = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";
// BLS12-381 scalar field modulus r
const char* bls12_381_scalar_modulus_hex = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

// Helper function to convert BigInt to byte array (little-endian)
void bigint_to_bytes_le(uint8_t *bytes, const BigInt &value, size_t len) {
    BigInt temp = value;
    for (size_t i = 0; i < len; i++) {
        bytes[i] = static_cast<uint8_t>(temp.to_ulong() & 0xff);
        temp = temp >> 8;
    }
}

// External BLST functions
extern "C" {
    extern const POINTonE1 BLS12_381_G1;
    extern const POINTonE2 BLS12_381_G2;
    void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
    void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
    void blst_p2_mult(POINTonE2 *out, const POINTonE2 *a, const byte *scalar, size_t nbits);
    void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a);
    void psi_wrapper_export(POINTonE2 *out, const POINTonE2 *in);
}

// Helper function to convert vec384 (Montgomery form) to BigInt
static BigInt vec384_montgomery_to_bigint_local(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    
    BigInt result(0);
    BigInt two_to_64 = BigInt(1);
    two_to_64 = two_to_64 << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(normal[i]));
    }
    return result;
}

// Helper function to convert POINTonE1_affine to our AffinePoint
template<class Ring>
AffinePoint<typename Ring::StandardElement> blst_g1_to_affine_point(
    const POINTonE1_affine &p,
    const Ring &ring
) {
    BigInt x = vec384_montgomery_to_bigint_local(p.X);
    BigInt y = vec384_montgomery_to_bigint_local(p.Y);
    
    AffinePoint<typename Ring::StandardElement> result;
    result.x = ring.from_bigint(x);
    result.y = ring.from_bigint(y);
    return result;
}

// Helper function to convert POINTonE2_affine to our AffinePoint
template<class Ring>
AffinePoint<typename FP2Ring<Ring>::StandardElement> blst_g2_to_affine_point(
    const POINTonE2_affine &p,
    const Ring &ring
) {
    FP2Ring<Ring> fp2ring(ring);
    auto x_pair = std::make_pair(
        vec384_montgomery_to_bigint_local(p.X[0]),
        vec384_montgomery_to_bigint_local(p.X[1])
    );
    auto y_pair = std::make_pair(
        vec384_montgomery_to_bigint_local(p.Y[0]),
        vec384_montgomery_to_bigint_local(p.Y[1])
    );
    
    AffinePoint<typename FP2Ring<Ring>::StandardElement> result;
    result.x = typename FP2Ring<Ring>::StandardElement(
        ring.from_bigint(x_pair.first),
        ring.from_bigint(x_pair.second)
    );
    result.y = typename FP2Ring<Ring>::StandardElement(
        ring.from_bigint(y_pair.first),
        ring.from_bigint(y_pair.second)
    );
    return result;
}

// Helper function to convert our ProjectivePoint to AffinePoint (BigInt)
template<class Ring>
AffinePoint<BigInt> to_bigint_g1(const ProjectivePoint<typename Ring::StandardElement> &point, const Ring &ring, const BigInt &modulus) {
    BigInt x = ring.to_bigint(point.x);
    BigInt y = ring.to_bigint(point.y);
    BigInt z = ring.to_bigint(point.z);
    
    // Check if point is at infinity (z == 0)
    if (z == BigInt(0)) {
        // Return point at infinity - in affine coordinates this is typically represented as (0, 0) or special handling
        // For comparison purposes, we'll use (0, 0) but this should not happen in normal scalar multiplication
        return AffinePoint<BigInt>{BigInt(0), BigInt(0)};
    }
    
    BigInt z_inv = z.mod_inverse(modulus);
    BigInt affine_x = (x * z_inv) % modulus;
    if (affine_x < 0) affine_x = affine_x + modulus;
    BigInt affine_y = (y * z_inv) % modulus;
    if (affine_y < 0) affine_y = affine_y + modulus;
    return AffinePoint<BigInt>{affine_x, affine_y};
}

// Helper function to convert our G2 ProjectivePoint to AffinePoint (BigInt pair)
template<class Ring>
AffinePoint<std::pair<BigInt, BigInt>> to_bigint_g2(const ProjectivePoint<typename FP2Ring<Ring>::StandardElement> &point, const Ring &ring, const BigInt &modulus) {
    FP2Ring<Ring> fp2ring(ring);
    BigInt x0 = ring.to_bigint(point.x.x());
    BigInt x1 = ring.to_bigint(point.x.y());
    BigInt y0 = ring.to_bigint(point.y.x());
    BigInt y1 = ring.to_bigint(point.y.y());
    BigInt z0 = ring.to_bigint(point.z.x());
    BigInt z1 = ring.to_bigint(point.z.y());
    
    // Check if point is at infinity (z == 0 + 0*i)
    if (z0 == BigInt(0) && z1 == BigInt(0)) {
        return AffinePoint<std::pair<BigInt, BigInt>>{
            std::make_pair(BigInt(0), BigInt(0)),
            std::make_pair(BigInt(0), BigInt(0))
        };
    }
    
    // Compute FP2 inverse: z = z0 + z1*i, z_inv = (z0 - z1*i) / (z0^2 + z1^2)
    // norm = z0^2 + z1^2
    BigInt norm = (z0 * z0 + z1 * z1) % modulus;
    if (norm < 0) norm = norm + modulus;
    BigInt norm_inv = norm.mod_inverse(modulus);
    
    // z_inv = (z0 * norm_inv, -z1 * norm_inv)
    BigInt z_inv0 = (z0 * norm_inv) % modulus;
    if (z_inv0 < 0) z_inv0 = z_inv0 + modulus;
    BigInt z_inv1 = (-z1 * norm_inv) % modulus;
    if (z_inv1 < 0) z_inv1 = z_inv1 + modulus;
    
    // Multiply x by z_inv: (x0 + x1*i) * (z_inv0 + z_inv1*i) = (x0*z_inv0 - x1*z_inv1) + (x0*z_inv1 + x1*z_inv0)*i
    BigInt affine_x0 = (x0 * z_inv0 - x1 * z_inv1) % modulus;
    if (affine_x0 < 0) affine_x0 = affine_x0 + modulus;
    BigInt affine_x1 = (x0 * z_inv1 + x1 * z_inv0) % modulus;
    if (affine_x1 < 0) affine_x1 = affine_x1 + modulus;
    
    // Multiply y by z_inv: (y0 + y1*i) * (z_inv0 + z_inv1*i) = (y0*z_inv0 - y1*z_inv1) + (y0*z_inv1 + y1*z_inv0)*i
    BigInt affine_y0 = (y0 * z_inv0 - y1 * z_inv1) % modulus;
    if (affine_y0 < 0) affine_y0 = affine_y0 + modulus;
    BigInt affine_y1 = (y0 * z_inv1 + y1 * z_inv0) % modulus;
    if (affine_y1 < 0) affine_y1 = affine_y1 + modulus;
    
    return AffinePoint<std::pair<BigInt, BigInt>>{
        std::make_pair(affine_x0, affine_x1),
        std::make_pair(affine_y0, affine_y1)
    };
}

// Helper function to compare G1 points
bool compare_g1_points(
    const AffinePoint<BigInt> &blst_result,
    const AffinePoint<BigInt> &our_result,
    const std::string &test_name
) {
    bool match = (blst_result.x == our_result.x && blst_result.y == our_result.y);
    
    if (!match) {
        std::cout << "✗ Mismatch in " << test_name << std::endl;
        std::cout << "  BLST X:  " << blst_result.x.to_string(16) << std::endl;
        std::cout << "  Ours X:  " << our_result.x.to_string(16) << std::endl;
        std::cout << "  BLST Y:  " << blst_result.y.to_string(16) << std::endl;
        std::cout << "  Ours Y:  " << our_result.y.to_string(16) << std::endl;
    }
    
    return match;
}

// Helper function to compare G2 points
bool compare_g2_points(
    const AffinePoint<std::pair<BigInt, BigInt>> &blst_result,
    const AffinePoint<std::pair<BigInt, BigInt>> &our_result,
    const std::string &test_name
) {
    bool match = (blst_result.x.first == our_result.x.first && 
                  blst_result.x.second == our_result.x.second &&
                  blst_result.y.first == our_result.y.first &&
                  blst_result.y.second == our_result.y.second);
    
    if (!match) {
        std::cout << "✗ Mismatch in " << test_name << std::endl;
        std::cout << "  BLST X:  (" << blst_result.x.first.to_string(16) << ", " << blst_result.x.second.to_string(16) << ")" << std::endl;
        std::cout << "  Ours X:  (" << our_result.x.first.to_string(16) << ", " << our_result.x.second.to_string(16) << ")" << std::endl;
        std::cout << "  BLST Y:  (" << blst_result.y.first.to_string(16) << ", " << blst_result.y.second.to_string(16) << ")" << std::endl;
        std::cout << "  Ours Y:  (" << our_result.y.first.to_string(16) << ", " << our_result.y.second.to_string(16) << ")" << std::endl;
    }
    
    return match;
}

// Helper function to convert our G2 ProjectivePoint to BLST POINTonE2 (Montgomery form)
template<class Ring>
void g2_projective_to_blst(const ProjectivePoint<typename FP2Ring<Ring>::StandardElement> &point, POINTonE2 &out, const Ring &ring) {
    // Convert to BigInt first
    BigInt x0 = ring.to_bigint(point.x.x());
    BigInt x1 = ring.to_bigint(point.x.y());
    BigInt y0 = ring.to_bigint(point.y.x());
    BigInt y1 = ring.to_bigint(point.y.y());
    BigInt z0 = ring.to_bigint(point.z.x());
    BigInt z1 = ring.to_bigint(point.z.y());
    
    // Ensure values are in range [0, q)
    BigInt q(bls12_381_modulus_hex, 16);
    x0 = x0 % q; if (x0 < 0) x0 = x0 + q;
    x1 = x1 % q; if (x1 < 0) x1 = x1 + q;
    y0 = y0 % q; if (y0 < 0) y0 = y0 + q;
    y1 = y1 % q; if (y1 < 0) y1 = y1 + q;
    z0 = z0 % q; if (z0 < 0) z0 = z0 + q;
    z1 = z1 % q; if (z1 < 0) z1 = z1 + q;
    
    // Convert to Montgomery form
    bigint_to_vec384_montgomery(out.X[0], x0);
    bigint_to_vec384_montgomery(out.X[1], x1);
    bigint_to_vec384_montgomery(out.Y[0], y0);
    bigint_to_vec384_montgomery(out.Y[1], y1);
    bigint_to_vec384_montgomery(out.Z[0], z0);
    bigint_to_vec384_montgomery(out.Z[1], z1);
}

// Helper function to convert BLST POINTonE2 to our G2 ProjectivePoint
template<class Ring>
ProjectivePoint<typename FP2Ring<Ring>::StandardElement> blst_to_g2_projective(const POINTonE2 &point, const Ring &ring) {
    FP2Ring<Ring> fp2ring(ring);
    
    // Convert from Montgomery form to BigInt
    BigInt x0 = vec384_montgomery_to_bigint_local(point.X[0]);
    BigInt x1 = vec384_montgomery_to_bigint_local(point.X[1]);
    BigInt y0 = vec384_montgomery_to_bigint_local(point.Y[0]);
    BigInt y1 = vec384_montgomery_to_bigint_local(point.Y[1]);
    BigInt z0 = vec384_montgomery_to_bigint_local(point.Z[0]);
    BigInt z1 = vec384_montgomery_to_bigint_local(point.Z[1]);
    
    // Convert to our format
    typename FP2Ring<Ring>::StandardElement x(
        ring.from_bigint(x0),
        ring.from_bigint(x1)
    );
    typename FP2Ring<Ring>::StandardElement y(
        ring.from_bigint(y0),
        ring.from_bigint(y1)
    );
    typename FP2Ring<Ring>::StandardElement z(
        ring.from_bigint(z0),
        ring.from_bigint(z1)
    );
    
    return ProjectivePoint<typename FP2Ring<Ring>::StandardElement>(x, y, z);
}

// Test GLV scalar multiplication for G1 (byte array version)
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_glv_g1_bytearray_config(const std::string &config_name) {
    std::cout << "=== Testing GLV Scalar Multiplication G1 (Byte Array) (" << config_name << "): BLST vs Our Implementation ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    using StandardElement = typename RingType::StandardElement;
    
    G1<RingType> g1_curve;
    GLVMult<RingType> glv_mult(g1_curve, ring);
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 10;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    // BLS12-381 scalar field modulus r
    BigInt r(bls12_381_scalar_modulus_hex, 16);
    
    std::cout << "Testing " << num_tests << " random scalar multiplications..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random scalar in range [0, r)
        BigInt scalar_bigint = BigInt::random(256) % r;
        
        // Convert to byte array (little-endian)
        byte scalar[32] = {0};
        bigint_to_bytes_le(scalar, scalar_bigint, 32);
        
        // Generate random G1 point using blst
        POINTonE1 random_point_proj;
        BigInt point_scalar_bigint = BigInt::random(256) % r;
        byte point_scalar[32] = {0};
        bigint_to_bytes_le(point_scalar, point_scalar_bigint, 32);
        blst_p1_mult(&random_point_proj, &BLS12_381_G1, point_scalar, 256);
        POINTonE1_affine random_point_affine;
        blst_p1_to_affine(&random_point_affine, &random_point_proj);
        
        // Compute scalar multiplication using BLST
        POINTonE1 blst_result_proj;
        blst_p1_mult(&blst_result_proj, &random_point_proj, scalar, 256);
        POINTonE1_affine blst_result_affine;
        blst_p1_to_affine(&blst_result_affine, &blst_result_proj);
        
        // Convert BLST result to BigInt
        AffinePoint<BigInt> blst_result_bigint{
            vec384_montgomery_to_bigint_local(blst_result_affine.X),
            vec384_montgomery_to_bigint_local(blst_result_affine.Y)
        };
        
        // Convert random point to our format
        auto P_affine = blst_g1_to_affine_point(random_point_affine, ring);
        typename GLVMult<RingType>::ProjectivePoint P_proj(P_affine.x, P_affine.y, ring.one());
        
        // Compute scalar multiplication using our GLV implementation (byte array version)
        std::array<uint8_t, 32> scalar_array;
        for (int i = 0; i < 32; i++) {
            scalar_array[i] = scalar[i];
        }
        typename GLVMult<RingType>::ProjectivePoint our_result_proj = glv_mult.multiply(P_proj, scalar_array, ring);
        
        // Convert our result to BigInt
        AffinePoint<BigInt> our_result_bigint = to_bigint_g1(our_result_proj, ring, q);
        
        // Compare results
        bool match = compare_g1_points(blst_result_bigint, our_result_bigint, 
                                       "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== GLV G1 (Byte Array) Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All GLV G1 (Byte Array) tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some GLV G1 (Byte Array) tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

// Test GLV scalar multiplication for G1 (BigInt version)
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_glv_g1_bigint_config(const std::string &config_name) {
    std::cout << "=== Testing GLV Scalar Multiplication G1 (BigInt) (" << config_name << "): BLST vs Our Implementation ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    using StandardElement = typename RingType::StandardElement;
    
    G1<RingType> g1_curve;
    GLVMult<RingType> glv_mult(g1_curve, ring);
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 10;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    // BLS12-381 scalar field modulus r
    BigInt r(bls12_381_scalar_modulus_hex, 16);
    
    std::cout << "Testing " << num_tests << " random scalar multiplications..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random scalar in range [0, r)
        BigInt scalar_bigint = BigInt::random(256) % r;
        
        // Convert to byte array (little-endian) for BLST
        byte scalar[32] = {0};
        bigint_to_bytes_le(scalar, scalar_bigint, 32);
        
        // Generate random G1 point using blst
        POINTonE1 random_point_proj;
        BigInt point_scalar_bigint = BigInt::random(256) % r;
        byte point_scalar[32] = {0};
        bigint_to_bytes_le(point_scalar, point_scalar_bigint, 32);
        blst_p1_mult(&random_point_proj, &BLS12_381_G1, point_scalar, 256);
        POINTonE1_affine random_point_affine;
        blst_p1_to_affine(&random_point_affine, &random_point_proj);
        
        // Compute scalar multiplication using BLST
        POINTonE1 blst_result_proj;
        blst_p1_mult(&blst_result_proj, &random_point_proj, scalar, 256);
        POINTonE1_affine blst_result_affine;
        blst_p1_to_affine(&blst_result_affine, &blst_result_proj);
        
        // Convert BLST result to BigInt
        AffinePoint<BigInt> blst_result_bigint{
            vec384_montgomery_to_bigint_local(blst_result_affine.X),
            vec384_montgomery_to_bigint_local(blst_result_affine.Y)
        };
        
        // Convert random point to our format
        auto P_affine = blst_g1_to_affine_point(random_point_affine, ring);
        typename GLVMult<RingType>::ProjectivePoint P_proj(P_affine.x, P_affine.y, ring.one());
        
        // Compute scalar multiplication using our GLV implementation (BigInt version)
        typename GLVMult<RingType>::ProjectivePoint our_result_proj = glv_mult.multiply(P_proj, scalar_bigint, ring);
        
        // Convert our result to BigInt
        AffinePoint<BigInt> our_result_bigint = to_bigint_g1(our_result_proj, ring, q);
        
        // Compare results
        bool match = compare_g1_points(blst_result_bigint, our_result_bigint, 
                                       "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== GLV G1 (BigInt) Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All GLV G1 (BigInt) tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some GLV G1 (BigInt) tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

// Test scalar_mult_table for G1
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_scalar_mult_table_g1_config(const std::string &config_name) {
    std::cout << "=== Testing Scalar Multiplication Table G1 (" << config_name << "): BLST vs Our Implementation ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    using StandardElement = typename RingType::StandardElement;
    
    G1<RingType> g1_curve;
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 10;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    std::cout << "Testing " << num_tests << " random scalar multiplications..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random 256-bit scalar
        byte scalar[32];
        for (int i = 0; i < 32; i++) {
            scalar[i] = static_cast<byte>(gen() & 0xff);
        }
        
        // Generate random G1 point using blst
        POINTonE1 random_point_proj;
        byte point_scalar[32];
        for (int i = 0; i < 32; i++) {
            point_scalar[i] = static_cast<byte>(gen() & 0xff);
        }
        blst_p1_mult(&random_point_proj, &BLS12_381_G1, point_scalar, 256);
        POINTonE1_affine random_point_affine;
        blst_p1_to_affine(&random_point_affine, &random_point_proj);
        
        // Compute scalar multiplication using BLST
        POINTonE1 blst_result_proj;
        blst_p1_mult(&blst_result_proj, &random_point_proj, scalar, 256);
        POINTonE1_affine blst_result_affine;
        blst_p1_to_affine(&blst_result_affine, &blst_result_proj);
        
        // Convert BLST result to BigInt
        AffinePoint<BigInt> blst_result_bigint{
            vec384_montgomery_to_bigint_local(blst_result_affine.X),
            vec384_montgomery_to_bigint_local(blst_result_affine.Y)
        };
        
        // Convert random point to our format
        auto P_affine = blst_g1_to_affine_point(random_point_affine, ring);
        typename G1<RingType>::ProjPoint P_proj(P_affine.x, P_affine.y, ring.one());
        
        // Convert scalar from bytes to BigInt (little-endian, matching BLST)
        BigInt scalar_bigint(0);
        for (int i = 31; i >= 0; i--) {
            scalar_bigint = scalar_bigint * BigInt(256) + BigInt(static_cast<unsigned char>(scalar[i]));
        }
        
        // Compute scalar multiplication using our scalar_mult_table implementation
        typename G1<RingType>::ProjPoint our_result_proj = scalar_mult_table<G1<RingType>, RingType, 256, 4>(
            scalar_bigint, P_proj, g1_curve, ring
        );
        
        // Convert our result to BigInt
        AffinePoint<BigInt> our_result_bigint = to_bigint_g1(our_result_proj, ring, q);
        
        // Compare results
        bool match = compare_g1_points(blst_result_bigint, our_result_bigint, 
                                       "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== Scalar Multiplication Table G1 Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All Scalar Multiplication Table G1 tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some Scalar Multiplication Table G1 tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

// Test scalar_mult_table with specific values from glv.py
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377, int scalar_bits=256>
void test_scalar_mult_table_glv_values(const std::string &config_name) {
    std::cout << "=== Testing Scalar Multiplication Table G1 with GLV.py values (" << config_name << ", " << scalar_bits << " bits) ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    
    G1<RingType> g1_curve;
    
    // BLS12-381 scalar field modulus r
    const char* bls12_381_scalar_modulus_hex = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
    BigInt r(bls12_381_scalar_modulus_hex, 16);
    
    // Test values from glv.py
    // Note: 2^256 (257 bits) is commented out as it's an edge case that BLST handles differently
    // (BLST returns point at infinity for this case, likely due to overflow handling)
    std::vector<BigInt> test_scalars = {
        BigInt(10),
        BigInt(28),
        BigInt(38),
        BigInt(64),
        r - BigInt(1),  // p - 1 (255 bits)
        (BigInt(1) << 256) - BigInt(1),  // 2^256 - 1 (256 bits)
        (BigInt(1) << 255) - BigInt(1),  // 2^255 - 1 (255 bits)
        // BigInt(1) << 256,  // 2^256 (257 bits) - commented out: BLST returns (0,0) for this case
        BigInt(0)
    };
    
    // Use generator G1 as the base point
    POINTonE1_affine g1_affine;
    blst_p1_to_affine(&g1_affine, &BLS12_381_G1);
    auto P_affine = blst_g1_to_affine_point(g1_affine, ring);
    typename G1<RingType>::ProjPoint P_proj(P_affine.x, P_affine.y, ring.one());
    
    int passed = 0;
    int failed = 0;
    
    // Filter test scalars based on bit length to match scalar_bits parameter
    std::vector<BigInt> filtered_scalars;
    for (size_t t = 0; t < test_scalars.size(); t++) {
        BigInt scalar = test_scalars[t];
        int scalar_bit_length = 0;
        BigInt temp_bits = scalar;
        while (temp_bits > 0) {
            scalar_bit_length++;
            temp_bits = temp_bits >> 1;
        }
        // Only include scalars that match the target bit length (within 1 bit tolerance)
        // For 255-bit test: include scalars <= 255 bits
        // For 256-bit test: include scalars <= 256 bits
        if (scalar_bits == 255 && scalar_bit_length <= 255) {
            filtered_scalars.push_back(scalar);
        } else if (scalar_bits == 256 && scalar_bit_length <= 256) {
            filtered_scalars.push_back(scalar);
        }
    }
    
    std::cout << "Testing " << filtered_scalars.size() << " specific scalar values from glv.py (filtered for " << scalar_bits << " bits)..." << std::endl << std::endl;
    
    for (size_t t = 0; t < filtered_scalars.size(); t++) {
        BigInt scalar = filtered_scalars[t];
        BigInt original_scalar = scalar;
        
        // For 256-bit scalars, don't reduce modulo r to match BLST behavior
        // BLST doesn't pre-reduce scalars that are >= r
        int scalar_bit_length = 0;
        BigInt temp_bits = scalar;
        while (temp_bits > 0) {
            scalar_bit_length++;
            temp_bits = temp_bits >> 1;
        }
        if (scalar_bit_length <= 255) {
            // Reduce scalar modulo r (scalar field modulus) for scalars <= 255 bits
            scalar = scalar % r;
        }
        // For 256+ bit scalars, use as-is (BLST will handle it)
        
        // Convert scalar to byte array (little-endian) for BLST
        byte scalar_bytes[32] = {0};
        BigInt temp = scalar;
        for (int i = 0; i < 32 && temp > 0; i++) {
            scalar_bytes[i] = static_cast<byte>(temp.to_ulong() & 0xff);
            temp = temp >> 8;
        }
        
        // Compute scalar multiplication using BLST
        POINTonE1 blst_result_proj;
        blst_p1_mult(&blst_result_proj, &BLS12_381_G1, scalar_bytes, scalar_bits);
        POINTonE1_affine blst_result_affine;
        blst_p1_to_affine(&blst_result_affine, &blst_result_proj);
        
        // Convert BLST result to BigInt
        AffinePoint<BigInt> blst_result_bigint{
            vec384_montgomery_to_bigint_local(blst_result_affine.X),
            vec384_montgomery_to_bigint_local(blst_result_affine.Y)
        };
        
        // Compute scalar multiplication using our scalar_mult_table implementation
        typename G1<RingType>::ProjPoint our_result_proj = scalar_mult_table<G1<RingType>, RingType, scalar_bits, 4>(
            scalar, P_proj, g1_curve, ring
        );
        
        // Convert our result to BigInt
        AffinePoint<BigInt> our_result_bigint = to_bigint_g1(our_result_proj, ring, q);
        
        // Compare results
        std::string scalar_str = scalar.to_string(16);
        bool match = compare_g1_points(blst_result_bigint, our_result_bigint, 
                                       "Scalar: 0x" + scalar_str);
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " (scalar: 0x" << scalar_str << ") passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " (scalar: 0x" << scalar_str << ") failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== GLV.py Values Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All GLV.py value tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some GLV.py value tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

// Test minus_psi against psi_wrapper_export
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_minus_psi(const std::string &config_name) {
    std::cout << "=== Testing minus_psi vs psi_wrapper_export (" << config_name << ") ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    using FP2RingType = FP2Ring<RingType>;
    
    FP2RingType fp2ring(ring);
    GLSMult<RingType> gls_mult(ring);
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 10;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    std::cout << "Testing " << num_tests << " random points..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random G2 point using blst
        POINTonE2 random_point_proj;
        byte point_scalar[32];
        for (int i = 0; i < 32; i++) {
            point_scalar[i] = static_cast<byte>(gen() & 0xff);
        }
        blst_p2_mult(&random_point_proj, &BLS12_381_G2, point_scalar, 256);
        
        // Convert to our format
        auto P_proj = blst_to_g2_projective(random_point_proj, ring);
        
        // Apply minus_psi (our implementation)
        typename GLSMult<RingType>::ProjectivePoint our_result = gls_mult.minus_psi(P_proj, fp2ring, ring);
        
        // Apply psi_wrapper_export (BLST)
        POINTonE2 blst_result;
        psi_wrapper_export(&blst_result, &random_point_proj);
        
        // Convert both results to BigInt for comparison
        AffinePoint<std::pair<BigInt, BigInt>> our_result_bigint = to_bigint_g2(our_result, ring, q);
        AffinePoint<std::pair<BigInt, BigInt>> blst_result_bigint = to_bigint_g2(blst_to_g2_projective(blst_result, ring), ring, q);
        
        // Compare results
        bool match = compare_g2_points(blst_result_bigint, our_result_bigint, 
                                       "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== minus_psi Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All minus_psi tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some minus_psi tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

// Test GLS scalar multiplication for G2 (byte array version)
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_gls_g2_bytearray_config(const std::string &config_name) {
    std::cout << "=== Testing GLS Scalar Multiplication G2 (Byte Array) (" << config_name << "): BLST vs Our Implementation ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    using FP2RingType = FP2Ring<RingType>;
    
    FP2RingType fp2ring(ring);
    GLSMult<RingType> gls_mult(ring);
    typename GLSMult<RingType>::G2Curve g2_curve;
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 10;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    // BLS12-381 scalar field modulus r
    BigInt r(bls12_381_scalar_modulus_hex, 16);
    
    std::cout << "Testing " << num_tests << " random scalar multiplications..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random scalar in range [0, r)
        BigInt scalar_bigint = BigInt::random(256) % r;
        
        // Convert to byte array (little-endian)
        byte scalar[32] = {0};
        bigint_to_bytes_le(scalar, scalar_bigint, 32);
        
        // Generate random G2 point using blst
        POINTonE2 random_point_proj;
        BigInt point_scalar_bigint = BigInt::random(256) % r;
        byte point_scalar[32] = {0};
        bigint_to_bytes_le(point_scalar, point_scalar_bigint, 32);
        blst_p2_mult(&random_point_proj, &BLS12_381_G2, point_scalar, 256);
        POINTonE2_affine random_point_affine;
        blst_p2_to_affine(&random_point_affine, &random_point_proj);
        
        // Compute scalar multiplication using BLST
        POINTonE2 blst_result_proj;
        blst_p2_mult(&blst_result_proj, &random_point_proj, scalar, 256);
        POINTonE2_affine blst_result_affine;
        blst_p2_to_affine(&blst_result_affine, &blst_result_proj);
        
        // Convert BLST result to BigInt
        AffinePoint<std::pair<BigInt, BigInt>> blst_result_bigint{
            std::make_pair(
                vec384_montgomery_to_bigint_local(blst_result_affine.X[0]),
                vec384_montgomery_to_bigint_local(blst_result_affine.X[1])
            ),
            std::make_pair(
                vec384_montgomery_to_bigint_local(blst_result_affine.Y[0]),
                vec384_montgomery_to_bigint_local(blst_result_affine.Y[1])
            )
        };
        
        // Convert random point to our format
        auto P_affine = blst_g2_to_affine_point(random_point_affine, ring);
        typename GLSMult<RingType>::ProjectivePoint P_proj(P_affine.x, P_affine.y, fp2ring.one());
        
        // Compute scalar multiplication using our GLS implementation (byte array version)
        std::array<uint8_t, 32> scalar_array;
        for (int i = 0; i < 32; i++) {
            scalar_array[i] = scalar[i];
        }
        typename GLSMult<RingType>::ProjectivePoint our_result_proj = gls_mult.multiply(P_proj, scalar_array, fp2ring, g2_curve, ring);
        
        // Convert our result to BigInt
        AffinePoint<std::pair<BigInt, BigInt>> our_result_bigint = to_bigint_g2(our_result_proj, ring, q);
        
        // Compare results
        bool match = compare_g2_points(blst_result_bigint, our_result_bigint, 
                                       "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== GLS G2 (Byte Array) Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All GLS G2 (Byte Array) tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some GLS G2 (Byte Array) tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

// Test GLS scalar multiplication for G2 (BigInt version)
template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_gls_g2_bigint_config(const std::string &config_name) {
    std::cout << "=== Testing GLS Scalar Multiplication G2 (BigInt) (" << config_name << "): BLST vs Our Implementation ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12> ring(q);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS, 12>;
    using FP2RingType = FP2Ring<RingType>;
    
    FP2RingType fp2ring(ring);
    GLSMult<RingType> gls_mult(ring);
    typename GLSMult<RingType>::G2Curve g2_curve;
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 10;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    // BLS12-381 scalar field modulus r
    BigInt r(bls12_381_scalar_modulus_hex, 16);
    
    std::cout << "Testing " << num_tests << " random scalar multiplications..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random scalar in range [0, r)
        BigInt scalar_bigint = BigInt::random(256) % r;
        
        // Convert to byte array (little-endian) for BLST
        byte scalar[32] = {0};
        bigint_to_bytes_le(scalar, scalar_bigint, 32);
        
        // Generate random G2 point using blst
        POINTonE2 random_point_proj;
        BigInt point_scalar_bigint = BigInt::random(256) % r;
        byte point_scalar[32] = {0};
        bigint_to_bytes_le(point_scalar, point_scalar_bigint, 32);
        blst_p2_mult(&random_point_proj, &BLS12_381_G2, point_scalar, 256);
        POINTonE2_affine random_point_affine;
        blst_p2_to_affine(&random_point_affine, &random_point_proj);
        
        // Compute scalar multiplication using BLST
        POINTonE2 blst_result_proj;
        blst_p2_mult(&blst_result_proj, &random_point_proj, scalar, 256);
        POINTonE2_affine blst_result_affine;
        blst_p2_to_affine(&blst_result_affine, &blst_result_proj);
        
        // Convert BLST result to BigInt
        AffinePoint<std::pair<BigInt, BigInt>> blst_result_bigint{
            std::make_pair(
                vec384_montgomery_to_bigint_local(blst_result_affine.X[0]),
                vec384_montgomery_to_bigint_local(blst_result_affine.X[1])
            ),
            std::make_pair(
                vec384_montgomery_to_bigint_local(blst_result_affine.Y[0]),
                vec384_montgomery_to_bigint_local(blst_result_affine.Y[1])
            )
        };
        
        // Convert random point to our format
        auto P_affine = blst_g2_to_affine_point(random_point_affine, ring);
        typename GLSMult<RingType>::ProjectivePoint P_proj(P_affine.x, P_affine.y, fp2ring.one());
        
        // Compute scalar multiplication using our GLS implementation (BigInt version)
        typename GLSMult<RingType>::ProjectivePoint our_result_proj = gls_mult.multiply(P_proj, scalar_bigint, fp2ring, g2_curve, ring);
        
        // Convert our result to BigInt
        AffinePoint<std::pair<BigInt, BigInt>> our_result_bigint = to_bigint_g2(our_result_proj, ring, q);
        
        // Compare results
        bool match = compare_g2_points(blst_result_bigint, our_result_bigint, 
                                       "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== GLS G2 (BigInt) Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All GLS G2 (BigInt) tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some GLS G2 (BigInt) tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    try {
        // GLV tests - test both byte array and BigInt versions
        // Test with 52-bit configuration
        test_glv_g1_bytearray_config<381, 8, 52, -1932, 2377>("52-bit");
        test_glv_g1_bigint_config<381, 8, 52, -1932, 2377>("52-bit");
        
        // Test with 50-bit configuration
        test_glv_g1_bytearray_config<381, 8, 50>("50-bit");
        test_glv_g1_bigint_config<381, 8, 50>("50-bit");

        /*
        
        // Test scalar_mult_table with GLV.py values (52-bit configuration, 255 bits)
        test_scalar_mult_table_glv_values<381, 8, 52, -1932, 2377, 255>("52-bit");
        
        // Test scalar_mult_table with GLV.py values (52-bit configuration, 256 bits)
        test_scalar_mult_table_glv_values<381, 8, 52, -1932, 2377, 256>("52-bit");
        
        // Test scalar_mult_table with GLV.py values (50-bit configuration, 255 bits)
        test_scalar_mult_table_glv_values<381, 8, 50, -1932, 2377, 255>("50-bit");
        
        // Test scalar_mult_table with GLV.py values (50-bit configuration, 256 bits)
        test_scalar_mult_table_glv_values<381, 8, 50, -1932, 2377, 256>("50-bit");
        
        // Test scalar_mult_table with random values (52-bit configuration)
        test_scalar_mult_table_g1_config<381, 8, 52, -1932, 2377>("52-bit");
        
        // Test scalar_mult_table with random values (50-bit configuration)
        test_scalar_mult_table_g1_config<381, 8, 50>("50-bit");
        */
        
        // Test minus_psi
        test_minus_psi<381, 8, 52, -1932, 2377>("52-bit");
        test_minus_psi<381, 8, 50>("50-bit");
        
        // GLS/G2 tests - test both byte array and BigInt versions
        test_gls_g2_bytearray_config<381, 8, 52, -1932, 2377>("52-bit");
        test_gls_g2_bigint_config<381, 8, 52, -1932, 2377>("52-bit");
        test_gls_g2_bytearray_config<381, 8, 50>("50-bit");
        test_gls_g2_bigint_config<381, 8, 50>("50-bit");
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
