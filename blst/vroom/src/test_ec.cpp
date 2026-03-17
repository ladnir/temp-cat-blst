#include "ec.hpp"
#include "bounded_ring.hpp"
#include "fp2.hpp"
#include "conversion_inversion.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <random>
#include <cassert>
extern "C" {
#include "../blst/vect.h"
#include "../blst/fields.h"
#include "../blst/consts.h"
#include "../blst/point.h"
}

// External BLST functions
extern "C" {
    extern const POINTonE1 BLS12_381_G1;
    extern const POINTonE2 BLS12_381_G2;
    void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
    void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
    void blst_p1_from_affine(POINTonE1 *out, const POINTonE1_affine *a);
    void blst_p1_add(POINTonE1 *out, const POINTonE1 *a, const POINTonE1 *b);
    void blst_p1_add_affine(POINTonE1 *out, const POINTonE1 *a, const POINTonE1_affine *b);
    void blst_p1_double(POINTonE1 *out, const POINTonE1 *a);
    void blst_p2_mult(POINTonE2 *out, const POINTonE2 *a, const byte *scalar, size_t nbits);
    void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a);
    void blst_p2_from_affine(POINTonE2 *out, const POINTonE2_affine *a);
    void blst_p2_add(POINTonE2 *out, const POINTonE2 *a, const POINTonE2 *b);
    void blst_p2_add_affine(POINTonE2 *out, const POINTonE2 *a, const POINTonE2_affine *b);
    void blst_p2_double(POINTonE2 *out, const POINTonE2 *a);
}

// Helper function to convert vec384x (Fp2, Montgomery form) to BigInt pair
static std::pair<BigInt, BigInt> vec384x_montgomery_to_bigint(const vec384x a) {
    return std::make_pair(
        vec384_montgomery_to_bigint(a[0]),
        vec384_montgomery_to_bigint(a[1])
    );
}

// Random generator
static std::random_device rd;
static std::mt19937_64 gen(rd());

void generate_random_point(POINTonE1_affine *point) {
    // Choose random 256 bit scalar
    byte scalar[32];
    for (int i = 0; i < 32; i++) {
        scalar[i] = static_cast<byte>(gen() & 0xff);
    }
    POINTonE1 p1_proj;
    blst_p1_mult(&p1_proj, &BLS12_381_G1, scalar, 256);
    blst_p1_to_affine(point, &p1_proj);
}

void generate_random_point(POINTonE2_affine *point) {
    // Choose random 256 bit scalar
    byte scalar[32];
    for (int i = 0; i < 32; i++) {
        scalar[i] = static_cast<byte>(gen() & 0xff);
    }
    POINTonE2 p2_proj;
    blst_p2_mult(&p2_proj, &BLS12_381_G2, scalar, 256);
    blst_p2_to_affine(point, &p2_proj);
}

template<class Ring>
AffinePoint<typename Ring::StandardElement> from_blst_g1(const POINTonE1_affine &point, const Ring &ring) {
    BigInt x = vec384_montgomery_to_bigint(point.X);
    BigInt y = vec384_montgomery_to_bigint(point.Y);
    return AffinePoint<typename Ring::StandardElement>{ring.from_bigint(x), ring.from_bigint(y)};
}

template<class Ring>
AffinePoint<typename FP2Ring<Ring>::StandardElement> from_blst_g2(const POINTonE2_affine &point, const Ring &ring) {
    FP2Ring<Ring> fp2ring(ring);
    auto x_pair = vec384x_montgomery_to_bigint(point.X);
    auto y_pair = vec384x_montgomery_to_bigint(point.Y);
    return AffinePoint<typename FP2Ring<Ring>::StandardElement>{
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(x_pair.first),
            ring.from_bigint(x_pair.second)
        ),
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(y_pair.first),
            ring.from_bigint(y_pair.second)
        )
    };
}

template<class Ring>
AffinePoint<BigInt> to_bigint_g1(const ProjectivePoint<typename Ring::StandardElement> &point, const Ring &ring, const BigInt &modulus) {
    BigInt x = ring.to_bigint(point.x);
    BigInt y = ring.to_bigint(point.y);
    BigInt z = ring.to_bigint(point.z);
    BigInt z_inv = z.mod_inverse(modulus);
    BigInt affine_x = (x * z_inv) % modulus;
    if (affine_x < 0) affine_x = affine_x + modulus;
    BigInt affine_y = (y * z_inv) % modulus;
    if (affine_y < 0) affine_y = affine_y + modulus;
    return AffinePoint<BigInt>{affine_x, affine_y};
}

template<class Ring>
AffinePoint<std::pair<BigInt, BigInt>> to_bigint_g2(const ProjectivePoint<typename FP2Ring<Ring>::StandardElement> &point, const Ring &ring, const BigInt &modulus) {
    FP2Ring<Ring> fp2ring(ring);
    BigInt x0 = ring.to_bigint(point.x.x());
    BigInt x1 = ring.to_bigint(point.x.y());
    BigInt y0 = ring.to_bigint(point.y.x());
    BigInt y1 = ring.to_bigint(point.y.y());
    BigInt z0 = ring.to_bigint(point.z.x());
    BigInt z1 = ring.to_bigint(point.z.y());
    
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

template<class Ring>
void test_g1_add(const Ring &ring, const G1<Ring> &curve, const BigInt &modulus) {
    POINTonE1_affine point1;
    POINTonE1_affine point2;
    generate_random_point(&point1);
    generate_random_point(&point2);

    POINTonE1 point1_proj;
    POINTonE1 point2_proj;
    blst_p1_from_affine(&point1_proj, &point1);
    blst_p1_from_affine(&point2_proj, &point2);
    POINTonE1 result_proj;
    blst_p1_add(&result_proj, &point1_proj, &point2_proj);
    POINTonE1_affine result_affine;
    blst_p1_to_affine(&result_affine, &result_proj);
   
    auto point1_affine = from_blst_g1(point1, ring);
    auto point2_affine = from_blst_g1(point2, ring);
    using ProjectivePoint = typename G1<Ring>::ProjPoint;
    ProjectivePoint point1_rns(point1_affine.x, point1_affine.y, ring.one());
    ProjectivePoint point2_rns(point2_affine.x, point2_affine.y, ring.one());
    ProjectivePoint result_rns = curve.add_point(point1_rns, point2_rns, ring);
    AffinePoint<BigInt> result_rns_bigint = to_bigint_g1(result_rns, ring, modulus);

    AffinePoint<BigInt> result_blst_bigint{vec384_montgomery_to_bigint(result_affine.X), vec384_montgomery_to_bigint(result_affine.Y)};
    assert(result_rns_bigint.x == result_blst_bigint.x && result_rns_bigint.y == result_blst_bigint.y);
}

template<class Ring>
void test_g1_mixed_add(const Ring &ring, const G1<Ring> &curve, const BigInt &modulus) {
    POINTonE1_affine point1_affine;
    POINTonE1_affine point2_affine;
    generate_random_point(&point1_affine);
    generate_random_point(&point2_affine);

    POINTonE1 point1_proj;
    blst_p1_from_affine(&point1_proj, &point1_affine);
    POINTonE1 result_proj;
    blst_p1_add_affine(&result_proj, &point1_proj, &point2_affine);
    POINTonE1_affine result_affine;
    blst_p1_to_affine(&result_affine, &result_proj);
   
    auto point1_aff = from_blst_g1(point1_affine, ring);
    auto point2_aff = from_blst_g1(point2_affine, ring);
    using ProjectivePoint = typename G1<Ring>::ProjPoint;
    ProjectivePoint point1_rns(point1_aff.x, point1_aff.y, ring.one());
    ProjectivePoint result_rns = curve.add_mixed_point(point1_rns, point2_aff, ring);
    AffinePoint<BigInt> result_rns_bigint = to_bigint_g1(result_rns, ring, modulus);

    AffinePoint<BigInt> result_blst_bigint{vec384_montgomery_to_bigint(result_affine.X), vec384_montgomery_to_bigint(result_affine.Y)};
    assert(result_rns_bigint.x == result_blst_bigint.x && result_rns_bigint.y == result_blst_bigint.y);
}

template<class Ring>
void test_g1_double(const Ring &ring, const G1<Ring> &curve, const BigInt &modulus) {
    POINTonE1_affine point_affine;
    generate_random_point(&point_affine);

    POINTonE1 point_proj;
    blst_p1_from_affine(&point_proj, &point_affine);
    POINTonE1 result_proj;
    blst_p1_double(&result_proj, &point_proj);
    POINTonE1_affine result_affine;
    blst_p1_to_affine(&result_affine, &result_proj);
   
    auto point_aff = from_blst_g1(point_affine, ring);
    using ProjectivePoint = typename G1<Ring>::ProjPoint;
    ProjectivePoint point_rns(point_aff.x, point_aff.y, ring.one());
    ProjectivePoint result_rns = curve.double_point(point_rns, ring);
    AffinePoint<BigInt> result_rns_bigint = to_bigint_g1(result_rns, ring, modulus);

    AffinePoint<BigInt> result_blst_bigint{vec384_montgomery_to_bigint(result_affine.X), vec384_montgomery_to_bigint(result_affine.Y)};
    assert(result_rns_bigint.x == result_blst_bigint.x && result_rns_bigint.y == result_blst_bigint.y);
}

template<class Ring>
void test_g2_add(const Ring &ring, const G2<Ring> &curve, const FP2Ring<Ring> &fp2ring, const BigInt &modulus) {
    POINTonE2_affine point1;
    POINTonE2_affine point2;
    generate_random_point(&point1);
    generate_random_point(&point2);

    POINTonE2 point1_proj;
    POINTonE2 point2_proj;
    blst_p2_from_affine(&point1_proj, &point1);
    blst_p2_from_affine(&point2_proj, &point2);
    POINTonE2 result_proj;
    blst_p2_add(&result_proj, &point1_proj, &point2_proj);
    POINTonE2_affine result_affine;
    blst_p2_to_affine(&result_affine, &result_proj);
   
    auto point1_aff = from_blst_g2(point1, ring);
    auto point2_aff = from_blst_g2(point2, ring);
    using ProjectivePoint = typename G2<Ring>::ProjPoint;
    ProjectivePoint point1_rns(point1_aff.x, point1_aff.y, fp2ring.one(ring));
    ProjectivePoint point2_rns(point2_aff.x, point2_aff.y, fp2ring.one(ring));
    ProjectivePoint result_rns = curve.add_point(point1_rns, point2_rns, fp2ring);
    AffinePoint<std::pair<BigInt, BigInt>> result_rns_bigint = to_bigint_g2(result_rns, ring, modulus);

    auto x_pair = vec384x_montgomery_to_bigint(result_affine.X);
    auto y_pair = vec384x_montgomery_to_bigint(result_affine.Y);
    AffinePoint<std::pair<BigInt, BigInt>> result_blst_bigint{
        std::make_pair(x_pair.first, x_pair.second),
        std::make_pair(y_pair.first, y_pair.second)
    };
    assert(result_rns_bigint.x.first == result_blst_bigint.x.first && result_rns_bigint.x.second == result_blst_bigint.x.second &&
           result_rns_bigint.y.first == result_blst_bigint.y.first && result_rns_bigint.y.second == result_blst_bigint.y.second);
}

template<class Ring>
void test_g2_mixed_add(const Ring &ring, const G2<Ring> &curve, const FP2Ring<Ring> &fp2ring, const BigInt &modulus) {
    POINTonE2_affine point1;
    POINTonE2_affine point2;
    generate_random_point(&point1);
    generate_random_point(&point2);

    POINTonE2 point1_proj;
    blst_p2_from_affine(&point1_proj, &point1);
    POINTonE2 result_proj;
    blst_p2_add_affine(&result_proj, &point1_proj, &point2);
    POINTonE2_affine result_affine;
    blst_p2_to_affine(&result_affine, &result_proj);
   
    auto point1_aff = from_blst_g2(point1, ring);
    auto point2_aff = from_blst_g2(point2, ring);
    using ProjectivePoint = typename G2<Ring>::ProjPoint;
    ProjectivePoint point1_rns(point1_aff.x, point1_aff.y, fp2ring.one(ring));
    ProjectivePoint result_rns = curve.add_mixed_point(point1_rns, point2_aff, fp2ring);
    AffinePoint<std::pair<BigInt, BigInt>> result_rns_bigint = to_bigint_g2(result_rns, ring, modulus);

    auto x_pair = vec384x_montgomery_to_bigint(result_affine.X);
    auto y_pair = vec384x_montgomery_to_bigint(result_affine.Y);
    AffinePoint<std::pair<BigInt, BigInt>> result_blst_bigint{
        std::make_pair(x_pair.first, x_pair.second),
        std::make_pair(y_pair.first, y_pair.second)
    };
    assert(result_rns_bigint.x.first == result_blst_bigint.x.first && result_rns_bigint.x.second == result_blst_bigint.x.second &&
           result_rns_bigint.y.first == result_blst_bigint.y.first && result_rns_bigint.y.second == result_blst_bigint.y.second);
}

template<class Ring>
void test_g2_double(const Ring &ring, const G2<Ring> &curve, const FP2Ring<Ring> &fp2ring, const BigInt &modulus) {
    POINTonE2_affine point;
    generate_random_point(&point);

    POINTonE2 point_proj;
    blst_p2_from_affine(&point_proj, &point);
    POINTonE2 result_proj;
    blst_p2_double(&result_proj, &point_proj);
    POINTonE2_affine result_affine;
    blst_p2_to_affine(&result_affine, &result_proj);
   
    auto point_aff = from_blst_g2(point, ring);
    using ProjectivePoint = typename G2<Ring>::ProjPoint;
    ProjectivePoint point_rns(point_aff.x, point_aff.y, fp2ring.one(ring));
    ProjectivePoint result_rns = curve.double_point(point_rns, fp2ring);
    AffinePoint<std::pair<BigInt, BigInt>> result_rns_bigint = to_bigint_g2(result_rns, ring, modulus);

    auto x_pair = vec384x_montgomery_to_bigint(result_affine.X);
    auto y_pair = vec384x_montgomery_to_bigint(result_affine.Y);
    AffinePoint<std::pair<BigInt, BigInt>> result_blst_bigint{
        std::make_pair(x_pair.first, x_pair.second),
        std::make_pair(y_pair.first, y_pair.second)
    };
    assert(result_rns_bigint.x.first == result_blst_bigint.x.first && result_rns_bigint.x.second == result_blst_bigint.x.second &&
           result_rns_bigint.y.first == result_blst_bigint.y.first && result_rns_bigint.y.second == result_blst_bigint.y.second);
}

template<class Ring>
void test_all(const Ring &ring, const BigInt &modulus) {
    FP2Ring<Ring> fp2ring(ring);
    G1<Ring> g1;
    G2<Ring> g2;
    for (int i = 0; i < 100; i++) {
        test_g1_add(ring, g1, modulus);
        test_g1_mixed_add(ring, g1, modulus);
        test_g1_double(ring, g1, modulus);
        test_g2_add(ring, g2, fp2ring, modulus);
        test_g2_mixed_add(ring, g2, fp2ring, modulus);
        test_g2_double(ring, g2, fp2ring, modulus);
    }
    std::cout << "All tests passed" << std::endl;
}

const char* bls12_381_modulus_hex = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";

int main() {
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<381, 8, 52, -1932, 2377, 12> bls_ring(q);
    test_all(bls_ring, q);
    BoundedRing<381, 8, 50, -1932, 2377, 12> bls_ring50(q);
    test_all(bls_ring50, q);
}