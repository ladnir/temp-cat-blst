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



	// Context for benchmarks with G1, G2 -> GT via pairing
	struct CoeffCtx
	{
		template<typename G, typename F>
		static auto plus(const G& a, const F& b) {
			G r;
			G::dadd(r, a, b);
			return r;
		}

		template<typename G, typename F>
		static auto minus(const G& a, const F& b) {
			G r;
			F neg_b;
			F::neg(neg_b, b);
			G::dadd(r, a, neg_b);
			return r;
		}

		// negation
		template<typename G>
		static auto neg(const G& a) {
			G r;
			G::neg(r, a);
			return r;
		}

		template<typename G, typename F>
		static auto mult(const G& a, const F& b) {
			if constexpr (std::is_same<G, F>::value) {
				G r;
				G::mul(r, a, b);
				return r;
			}
			else {
				G r;
				F b_non_mont;
				F::from_mont(b_non_mont, b);
				G::mult_vartime(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
				return r;
			}
		}

		static auto mult(const G1Point& a, const G2Point& b) {
			G1AffinePoint p_aff;
			G1Point::to_affine(p_aff, a);
			G2AffinePoint q_aff;
			G2Point::to_affine(q_aff, b);

			MillerResult miller;
			GT::miller_loop(miller, &q_aff, &p_aff);
			GT result;
			GT::final_exp(result, miller);
			return result;
		}

		static auto inverse(const Fr& a) {
			Fr r;
			Fr::reciprocal(r, a);
			return r;
		}

		template<typename F>
		static auto from_u64(uint64_t a) {
			Fr r;
			Fr::from_uint64(r, a);
			return r;
		}
	};
	using G1CoeffCtx = CoeffCtx;
	using G2CoeffCtx = CoeffCtx;
	using GTCoeffCtx = CoeffCtx;
}
