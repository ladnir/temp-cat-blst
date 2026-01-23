#pragma once
#include "blst/fields/Fr.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"

namespace blst
{





    struct FrCoeffCtx
    {
        static auto plus(const Fr& a, const Fr& b) {
            Fr r;
            Fr::add(r, a, b);
            return r;
        }
        static auto minus(const Fr& a, const Fr& b) {
            Fr r;
            Fr::sub(r, a, b);
            return r;
        }

        static auto mult(const Fr& a, const Fr& b) {
            Fr r;
            Fr::mul(r, a, b);
            return r;
        }
        static auto inverse(const Fr& a) {
            Fr r;
            Fr::reciprocal(r, a);
            return r;
        }

        template<typename F>
        static auto from_u64(u64 a) {
			return Fr::from_uint64(a);
        }

    };



    // Context for FFT operations on G1 points with Fr scalars (FFT in the exponent)
    struct G1CoeffCtx
    {
        // G1 + G1 -> G1 (point addition)
        static auto plus(const G1Point& a, const G1Point& b) {
            G1Point r;
            G1Point::dadd(r, a, b);
            return r;
        }

        // G1 - G1 -> G1 (point subtraction)
        static auto minus(const G1Point& a, const G1Point& b) {
            G1Point neg_b;
            G1Point::neg(neg_b, b);
            G1Point r;
            G1Point::dadd(r, a, neg_b);
            return r;
        }

        // G1 * Fr -> G1 (scalar multiplication)
        static auto mult(const G1Point& a, const Fr& b) {
            G1Point r;
            // Convert Fr from Montgomery form to scalar bytes for mult
            Fr b_non_mont;
            Fr::from_mont(b_non_mont, b);
            G1Point::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
            return r;
        }

        // Inverse of Fr scalar (for iFFT scaling)
        static auto inverse(const Fr& a) {
            Fr r;
            Fr::reciprocal(r, a);
            return r;
        }

        template<typename F>
        static auto from_u64(u64 a) {
            return Fr::from_uint64(a);
        }
    };



    // Context for FFT operations on G2 points with Fr scalars (FFT in the exponent)
    struct G2CoeffCtx
    {
        // G2 + G2 -> G2 (point addition)
        static auto plus(const G2Point& a, const G2Point& b) {
            G2Point r;
            G2Point::dadd(r, a, b);
            return r;
        }

        // G2 - G2 -> G2 (point subtraction)
        static auto minus(const G2Point& a, const G2Point& b) {
            G2Point neg_b;
            G2Point::neg(neg_b, b);
            G2Point r;
            G2Point::dadd(r, a, neg_b);
            return r;
        }

        // G2 * Fr -> G2 (scalar multiplication)
        static auto mult(const G2Point& a, const Fr& b) {
            G2Point r;
            // Convert Fr from Montgomery form to scalar bytes for mult
            Fr b_non_mont;
            Fr::from_mont(b_non_mont, b);
            G2Point::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
            return r;
        }

        // Inverse of Fr scalar (for iFFT scaling)
        static auto inverse(const Fr& a) {
            Fr r;
            Fr::reciprocal(r, a);
            return r;
        }

        template<typename F>
        static auto from_u64(u64 a) {
            return Fr::from_uint64(a);
        }
    };




    // Context for FFT operations on GT elements with Fr scalars
    struct GTCoeffCtx
    {
        // GT + GT -> GT (group addition = Fp12 multiplication)
        static auto plus(const GT& a, const GT& b) {
            GT r;
            GT::add(r, a, b);
            return r;
        }

        // GT - GT -> GT (group subtraction = Fp12 division)
        static auto minus(const GT& a, const GT& b) {
            GT r;
            GT::sub(r, a, b);
            return r;
        }

        // GT * Fr -> GT (scalar multiplication = Fp12 exponentiation)
        static auto mult(const GT& a, const Fr& b) {
            GT r;
            Fr b_non_mont;
            Fr::from_mont(b_non_mont, b);
            GT::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
            return r;
        }


        // Fr * Fr -> GT (scalar multiplication = Fp12 exponentiation)
        static auto mult(const Fr& a, const Fr& b) {
            Fr r;
            Fr::mul(r, a, b);
            return r;
        }

        // Inverse of Fr scalar (for iFFT scaling)
        static auto inverse(const Fr& a) {
            Fr r;
            Fr::reciprocal(r, a);
            return r;
        }

        template<typename F>
        static auto from_u64(u64 a) {
            return Fr::from_uint64(a);
        }
    };
}