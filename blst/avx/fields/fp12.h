// SPDX-License-Identifier: MIT
// Copyright 2026 Simon Langowski
// Copyright 2026 Category Labs
//
// Adapted from the VROOM Fp12 arithmetic at commit d6ecd534. CatBLST retains
// only the selected MatrixNoK schedules and makes their representation,
// reduction, timing, and alias contracts explicit.
#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include "blst/avx/fields/fp2.h"


#if defined(_MSC_VER)
#define BLST_AVX_FP12_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define BLST_AVX_FP12_NOINLINE __attribute__((noinline))
#else
#define BLST_AVX_FP12_NOINLINE
#endif

namespace blst::avx {

template<class Ring, class Constants>
class FP12Ring {

	public:
	using FP2StandardElement = typename FP2Ring<Ring>::StandardElement;

	// Fp12 is stored as six Fp2 coefficients in the scalar tower's established
	// order, flattened to twelve consecutive base-field coordinates. Bounds are
	// compile-time producer claims; constructing a CommonElement does not
	// validate either its field value or its physical/RNS interval.
	//
	// General multiply, square, Frobenius, conjugation, and sparse Miller
	// operations accept any valid Fp12 value. Cyclotomic square requires a
	// cyclotomic/GT input. Inversion requires nonzero input and inherits the
	// supplied inverter's timing; the selected inverter is explicitly variable
	// time. Operations are fixed-shape and allocation-free apart from that
	// provider contract, and perform no Release validation. Returned-value
	// operations naturally support exact input/output storage reuse by callers;
	// the one output-reference sparse multiply reads all inputs before assigning.
	template<int64_t common_bound>
	using CommonElement = std::array<typename Ring::template ReducedElement<common_bound>, 12>;

	using StandardElement = CommonElement<1>;
	
	// Shorthand for integral constants
	static constexpr auto _2 = std::integral_constant<int, 2>{};
	static constexpr auto _3 = std::integral_constant<int, 3>{};
	static constexpr auto _6 = std::integral_constant<int, 6>{};
	// The source-level bound algebra is conservative by one relative to the
	// proved correlated-value maxima (28 and 101 respectively).
	static constexpr int64_t first_cyclotomic_reduce_bound = 29;
	static constexpr int64_t second_cyclotomic_reduce_bound = 102;

	FP12Ring() noexcept = default;

	template<int64_t upper_bound, class E1, class E2, class E3, class E4, class E5, class E6>
	BLST_FORCEINLINE auto to_common(const Ring &ring, const E1 &e1, const E2 &e2, const E3 &e3, const E4 &e4, const E5 &e5, const E6 &e6) const {
		using ResultType = typename Ring::template SmallReducedElement<upper_bound>;
		return std::array<ResultType, 6>{{
			ring.prep_expand(e1).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e2).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e3).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e4).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e5).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e6).template recast_rns_bounds<0, upper_bound>(),
		}};
	}

	template<int64_t Z>
	BLST_FORCEINLINE auto recast_common(const CommonElement<Z> &x) const {
		constexpr int N = 6;
		return CommonElement<N> {
			x[0].template recast_rns_bounds<0, N>(),
			x[1].template recast_rns_bounds<0, N>(),
			x[2].template recast_rns_bounds<0, N>(),
			x[3].template recast_rns_bounds<0, N>(),
			x[4].template recast_rns_bounds<0, N>(),
			x[5].template recast_rns_bounds<0, N>(),
			x[6].template recast_rns_bounds<0, N>(),
			x[7].template recast_rns_bounds<0, N>(),
			x[8].template recast_rns_bounds<0, N>(),
			x[9].template recast_rns_bounds<0, N>(),
			x[10].template recast_rns_bounds<0, N>(),
			x[11].template recast_rns_bounds<0, N>(),
		};
	}

	template<int64_t Z>
	BLST_FORCEINLINE StandardElement narrow_to_standard(const CommonElement<Z> &x, const Ring &ring) const {
		const auto one = ring.one();
		return ring.batch_reduce_expand(
			ring.prep_left(x[0]) * one, ring.prep_left(x[1]) * one,
			ring.prep_left(x[2]) * one, ring.prep_left(x[3]) * one,
			ring.prep_left(x[4]) * one, ring.prep_left(x[5]) * one,
			ring.prep_left(x[6]) * one, ring.prep_left(x[7]) * one,
			ring.prep_left(x[8]) * one, ring.prep_left(x[9]) * one,
			ring.prep_left(x[10]) * one, ring.prep_left(x[11]) * one);
	}

	template<int64_t max_bound, typename... Args>
	BLST_FORCEINLINE auto batch_reduce(const Ring &ring, const Args&... args) const {
		constexpr int batch = static_cast<int>(sizeof...(Args));
		static_assert(batch > 0);
		const std::array<typename Ring::template ReadyElement<max_bound>, batch>
			ready{ring.template ready<max_bound>(args)...};
		return ring.template batch_reduce<batch, max_bound>(ready);
	}

	StandardElement one(const Ring &ring) const {
		return std::array<typename Ring::StandardElement, 12> {ring.one(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero()};
	}


	/* 
	Polynomial 0(22) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3 +1 x6y10 -1 x6y11 -1 x7y10 -1 x7y11 +1 x8y8 -1 x8y9 -1 x9y8 -1 x9y9 +1 x10y6 -1 x10y7 -1 x11y6 -1 x11y7

	Polynomial 1(22) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3 +1 x6y10 +1 x6y11 +1 x7y10 -1 x7y11 +1 x8y8 +1 x8y9 +1 x9y8 -1 x9y9 +1 x10y6 +1 x10y7 +1 x11y6 -1 x11y7

	Polynomial 2(18) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5 +1 x6y6 -1 x7y7 +1 x8y10 -1 x8y11 -1 x9y10 -1 x9y11 +1 x10y8 -1 x10y9 -1 x11y8 -1 x11y9

	Polynomial 3(18) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5 +1 x6y7 +1 x7y6 +1 x8y10 +1 x8y11 +1 x9y10 -1 x9y11 +1 x10y8 +1 x10y9 +1 x11y8 -1 x11y9

	Polynomial 4(14) : 1 x0y4 -1 x1y5 +1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1 +1 x6y8 -1 x7y9 +1 x8y6 -1 x9y7 +1 x10y10 -1 x10y11 -1 x11y10 -1 x11y11

	Polynomial 5(14) : 1 x0y5 +1 x1y4 +1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0 +1 x6y9 +1 x7y8 +1 x8y7 +1 x9y6 +1 x10y10 +1 x10y11 +1 x11y10 -1 x11y11

	Polynomial 6(20) : 1 x0y6 -1 x1y7 +1 x2y10 -1 x2y11 -1 x3y10 -1 x3y11 +1 x4y8 -1 x4y9 -1 x5y8 -1 x5y9 +1 x6y0 -1 x7y1 +1 x8y4 -1 x8y5 -1 x9y4 -1 x9y5 +1 x10y2 -1 x10y3 -1 x11y2 -1 x11y3

	Polynomial 7(20) : 1 x0y7 +1 x1y6 +1 x2y10 +1 x2y11 +1 x3y10 -1 x3y11 +1 x4y8 +1 x4y9 +1 x5y8 -1 x5y9 +1 x6y1 +1 x7y0 +1 x8y4 +1 x8y5 +1 x9y4 -1 x9y5 +1 x10y2 +1 x10y3 +1 x11y2 -1 x11y3

	Polynomial 8(16) : 1 x0y8 -1 x1y9 +1 x2y6 -1 x3y7 +1 x4y10 -1 x4y11 -1 x5y10 -1 x5y11 +1 x6y2 -1 x7y3 +1 x8y0 -1 x9y1 +1 x10y4 -1 x10y5 -1 x11y4 -1 x11y5

	Polynomial 9(16) : 1 x0y9 +1 x1y8 +1 x2y7 +1 x3y6 +1 x4y10 +1 x4y11 +1 x5y10 -1 x5y11 +1 x6y3 +1 x7y2 +1 x8y1 +1 x9y0 +1 x10y4 +1 x10y5 +1 x11y4 -1 x11y5

	Polynomial 10(12) : 1 x0y10 -1 x1y11 +1 x2y8 -1 x3y9 +1 x4y6 -1 x5y7 +1 x6y4 -1 x7y5 +1 x8y2 -1 x9y3 +1 x10y0 -1 x11y1

	Polynomial 11(12) : 1 x0y11 +1 x1y10 +1 x2y9 +1 x3y8 +1 x4y7 +1 x5y6 +1 x6y5 +1 x7y4 +1 x8y3 +1 x9y2 +1 x10y1 +1 x11y0
	*/

	private:
	// Compute one six-coordinate half of a dense Fp12 product at a time. The
	// split is algebraically neutral: together the helpers evaluate the twelve
	// polynomials above. Keeping the two reduction batches in separate
	// non-inlined frames bounds the live AVX-512 state and prevents the compiler
	// from spilling the monolithic polynomial schedule excessively.
	//
	// Preconditions: x and y contain valid Fp12 values satisfying their stated
	// compile-time bounds. result is fresh storage and does not alias either
	// input. Each helper writes only its named half of result. The schedule is
	// fixed-shape, allocation-free, and does not inspect field values.
	template<int64_t N1, int64_t N2>
	BLST_AVX_FP12_NOINLINE void fp12_flat_mul_first_split(
		StandardElement& result,
		const CommonElement<N1>& x,
		const CommonElement<N2>& y,
		const Ring& ring) const {
		const auto offset = Constants::offset();
		auto minus_y1 = ring.negate(y[1]);
		auto minus_y3 = ring.negate(y[3]);
		auto minus_y5 = ring.negate(y[5]);
		auto minus_y7 = ring.negate(y[7]);
		auto minus_y9 = ring.negate(y[9]);
		auto minus_y11 = ring.negate(y[11]);

		auto poly_0_pos_shared = offset + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3 + x[6] * y[10] + x[7] * minus_y11 + x[8] * y[8] + x[9] * minus_y9 + x[10] * y[6] + x[11] * minus_y7;
		auto poly_0_neg_shared = x[2] * y[5] + x[3] * y[4] + x[4] * y[3] + x[5] * y[2] + x[6] * y[11] + x[7] * y[10] + x[8] * y[9] + x[9] * y[8] + x[10] * y[7] + x[11] * y[6];
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
		auto poly_2_pos_shared = offset + x[4] * y[4] + x[5] * minus_y5 + x[8] * y[10] + x[9] * minus_y11 + x[10] * y[8] + x[11] * minus_y9;
		auto poly_2_neg_shared = x[4] * y[5] + x[5] * y[4] + x[8] * y[11] + x[9] * y[10] + x[10] * y[9] + x[11] * y[8];
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y[2] + x[1] * minus_y3 + x[2] * y[0] + x[3] * minus_y1 + x[6] * y[6] + x[7] * minus_y7;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0] + x[6] * y[7] + x[7] * y[6];
		auto poly_4_pos_shared = offset + x[10] * y[10] + x[11] * minus_y11;
		auto poly_4_neg_shared = x[10] * y[11] + x[11] * y[10];
		auto poly_4_unique = poly_4_pos_shared - poly_4_neg_shared + x[0] * y[4] + x[1] * minus_y5 + x[2] * y[2] + x[3] * minus_y3 + x[4] * y[0] + x[5] * minus_y1 + x[6] * y[8] + x[7] * minus_y9 + x[8] * y[6] + x[9] * minus_y7;
		auto poly_5_unique = poly_4_pos_shared + poly_4_neg_shared + x[0] * y[5] + x[1] * y[4] + x[2] * y[3] + x[3] * y[2] + x[4] * y[1] + x[5] * y[0] + x[6] * y[9] + x[7] * y[8] + x[8] * y[7] + x[9] * y[6];

		const auto expanded = ring.batch_reduce_expand(
			poly_0_unique, poly_1_unique, poly_2_unique,
			poly_3_unique, poly_4_unique, poly_5_unique);
		result[0] = expanded[0];
		result[1] = expanded[1];
		result[2] = expanded[2];
		result[3] = expanded[3];
		result[4] = expanded[4];
		result[5] = expanded[5];
	}

	template<int64_t N1, int64_t N2>
	BLST_AVX_FP12_NOINLINE void fp12_flat_mul_second_split(
		StandardElement& result,
		const CommonElement<N1>& x,
		const CommonElement<N2>& y,
		const Ring& ring) const {
		const auto offset = Constants::offset();
		auto minus_y1 = ring.negate(y[1]);
		auto minus_y3 = ring.negate(y[3]);
		auto minus_y5 = ring.negate(y[5]);
		auto minus_y7 = ring.negate(y[7]);
		auto minus_y9 = ring.negate(y[9]);
		auto minus_y11 = ring.negate(y[11]);

		auto poly_6_pos_shared = offset + x[2] * y[10] + x[3] * minus_y11 + x[4] * y[8] + x[5] * minus_y9 + x[8] * y[4] + x[9] * minus_y5 + x[10] * y[2] + x[11] * minus_y3;
		auto poly_6_neg_shared = x[2] * y[11] + x[3] * y[10] + x[4] * y[9] + x[5] * y[8] + x[8] * y[5] + x[9] * y[4] + x[10] * y[3] + x[11] * y[2];
		auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[0] * y[6] + x[1] * minus_y7 + x[6] * y[0] + x[7] * minus_y1;
		auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[0] * y[7] + x[1] * y[6] + x[6] * y[1] + x[7] * y[0];
		auto poly_8_pos_shared = offset + x[4] * y[10] + x[5] * minus_y11 + x[10] * y[4] + x[11] * minus_y5;
		auto poly_8_neg_shared = x[4] * y[11] + x[5] * y[10] + x[10] * y[5] + x[11] * y[4];
		auto poly_8_unique = poly_8_pos_shared - poly_8_neg_shared + x[0] * y[8] + x[1] * minus_y9 + x[2] * y[6] + x[3] * minus_y7 + x[6] * y[2] + x[7] * minus_y3 + x[8] * y[0] + x[9] * minus_y1;
		auto poly_9_unique = poly_8_pos_shared + poly_8_neg_shared + x[0] * y[9] + x[1] * y[8] + x[2] * y[7] + x[3] * y[6] + x[6] * y[3] + x[7] * y[2] + x[8] * y[1] + x[9] * y[0];
		auto poly_10_unique = x[0] * y[10] + x[1] * minus_y11 + x[2] * y[8] + x[3] * minus_y9 + x[4] * y[6] + x[5] * minus_y7 + x[6] * y[4] + x[7] * minus_y5 + x[8] * y[2] + x[9] * minus_y3 + x[10] * y[0] + x[11] * minus_y1;
		auto poly_11_unique = x[0] * y[11] + x[1] * y[10] + x[2] * y[9] + x[3] * y[8] + x[4] * y[7] + x[5] * y[6] + x[6] * y[5] + x[7] * y[4] + x[8] * y[3] + x[9] * y[2] + x[10] * y[1] + x[11] * y[0];

		const auto expanded = ring.batch_reduce_expand(
			poly_6_unique, poly_7_unique, poly_8_unique,
			poly_9_unique, poly_10_unique, poly_11_unique);
		result[6] = expanded[0];
		result[7] = expanded[1];
		result[8] = expanded[2];
		result[9] = expanded[3];
		result[10] = expanded[4];
		result[11] = expanded[5];
	}

	public:
	// Return x * y in Fp12. Inputs must be valid CommonElement values satisfying
	// their declared bounds. The returned value is fully reduced and expanded;
	// the operation is fixed-shape and allocation-free.
	template<int64_t N1, int64_t N2>
	inline StandardElement fp12_flat_mul(
		const CommonElement<N1>& x,
		const CommonElement<N2>& y,
		const Ring& ring) const {
		StandardElement result;
		fp12_flat_mul_first_split(result, x, y, ring);
		fp12_flat_mul_second_split(result, x, y, ring);
		return result;
	}

	/* 
	Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3 +1 x8y4 -1 x8y5 -1 x9y4 -1 x9y5

	Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3 +1 x8y4 +1 x8y5 +1 x9y4 -1 x9y5

	Polynomial 2(8) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x10y4 -1 x10y5 -1 x11y4 -1 x11y5

	Polynomial 3(8) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x10y4 +1 x10y5 +1 x11y4 -1 x11y5

	Polynomial 4(6) : 1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1 +1 x6y4 -1 x7y5

	Polynomial 5(6) : 1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0 +1 x6y5 +1 x7y4

	Polynomial 6(10) : 1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5 +1 x6y0 -1 x7y1 +1 x10y2 -1 x10y3 -1 x11y2 -1 x11y3

	Polynomial 7(10) : 1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5 +1 x6y1 +1 x7y0 +1 x10y2 +1 x10y3 +1 x11y2 -1 x11y3

	Polynomial 8(6) : 1 x0y4 -1 x1y5 +1 x6y2 -1 x7y3 +1 x8y0 -1 x9y1

	Polynomial 9(6) : 1 x0y5 +1 x1y4 +1 x6y3 +1 x7y2 +1 x8y1 +1 x9y0

	Polynomial 10(6) : 1 x2y4 -1 x3y5 +1 x8y2 -1 x9y3 +1 x10y0 -1 x11y1

	Polynomial 11(6) : 1 x2y5 +1 x3y4 +1 x8y3 +1 x9y2 +1 x10y1 +1 x11y0
	*/
	
	template<class CX, class CY>
        BLST_FORCEINLINE void flat_mulxy00z0(StandardElement &result, const StandardElement &x, const FP2<CX, CY> &constant_coord, const FP2StandardElement &x_coord, const FP2StandardElement &y_coord, const Ring &ring) const {
		auto y0 = constant_coord.x();
		auto y1 = constant_coord.y();
		auto y2 = x_coord.x();
		auto y3 = x_coord.y();
		auto y8 = y_coord.x();
		auto y9 = y_coord.y();
		auto minus_y1 = ring.negate(y1);
		auto minus_y3 = ring.negate(y3);
		auto minus_y9 = ring.negate(y9);
		// 0 1
		auto poly_0_pos_shared = Constants::offset() + x[4] * y2 + x[5] * minus_y3 + x[8] * y8 + x[9] * minus_y9;
		auto poly_0_neg_shared = x[4] * y3 + x[5] * y2 + x[8] * y9 + x[9] * y8;
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y0 + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y1 + x[1] * y0;
		// 2 3
		auto poly_2_pos_shared = Constants::offset() + x[10] * y8 + x[11] * minus_y9;
		auto poly_2_neg_shared = x[10] * y9 + x[11] * y8;
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y2 + x[1] * minus_y3 + x[2] * y0 + x[3] * minus_y1;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y3 + x[1] * y2 + x[2] * y1 + x[3] * y0;
		// 4 5
		auto poly_4_unique = x[2] * y2 + x[3] * minus_y3 + x[4] * y0 + x[5] * minus_y1 + x[6] * y8 + x[7] * minus_y9;
		auto poly_5_unique = x[2] * y3 + x[3] * y2 + x[4] * y1 + x[5] * y0 + x[6] * y9 + x[7] * y8;
		// 6 7
		auto poly_6_pos_shared = Constants::offset() + x[4] * y8 + x[5] * minus_y9 + x[10] * y2 + x[11] * minus_y3;
		auto poly_6_neg_shared = x[4] * y9 + x[5] * y8 + x[10] * y3 + x[11] * y2;
		auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[6] * y0 + x[7] * minus_y1;
		auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[6] * y1 + x[7] * y0;
		// 8 9
		auto poly_8_unique = x[0] * y8 + x[1] * minus_y9 + x[6] * y2 + x[7] * minus_y3 + x[8] * y0 + x[9] * minus_y1;
		auto poly_9_unique = x[0] * y9 + x[1] * y8 + x[6] * y3 + x[7] * y2 + x[8] * y1 + x[9] * y0;
		// 10 11
		auto poly_10_unique = x[2] * y8 + x[3] * minus_y9 + x[8] * y2 + x[9] * minus_y3 + x[10] * y0 + x[11] * minus_y1;
		auto poly_11_unique = x[2] * y9 + x[3] * y8 + x[8] * y3 + x[9] * y2 + x[10] * y1 + x[11] * y0;
		result = ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
	}
	/* 
	Polynomial 0(13) : 1 x0x0 -1 x1x1 +1 x8x8 -1 x9x9 +2 x2x4 -2 x2x5 -2 x3x4 -2 x3x5 +2 x6x10 -2 x6x11 -2 x7x10 -2 x7x11 -2 x8x9

	Polynomial 1(12) : 1 x8x8 -1 x9x9 +2 x0x1 +2 x2x4 +2 x2x5 +2 x3x4 -2 x3x5 +2 x6x10 +2 x6x11 +2 x7x10 -2 x7x11 +2 x8x9

	Polynomial 2(11) : 1 x4x4 -1 x5x5 +1 x6x6 -1 x7x7 +2 x0x2 -2 x1x3 -2 x4x5 +2 x8x10 -2 x8x11 -2 x9x10 -2 x9x11

	Polynomial 3(10) : 1 x4x4 -1 x5x5 +2 x0x3 +2 x1x2 +2 x4x5 +2 x6x7 +2 x8x10 +2 x8x11 +2 x9x10 -2 x9x11

	Polynomial 4(9) : 1 x2x2 -1 x3x3 +1 x10x10 -1 x11x11 +2 x0x4 -2 x1x5 +2 x6x8 -2 x7x9 -2 x10x11

	Polynomial 5(8) : 1 x10x10 -1 x11x11 +2 x0x5 +2 x1x4 +2 x2x3 +2 x6x9 +2 x7x8 +2 x10x11

	Polynomial 6(10) : 2 x0x6 -2 x1x7 +2 x2x10 -2 x2x11 -2 x3x10 -2 x3x11 +2 x4x8 -2 x4x9 -2 x5x8 -2 x5x9

	Polynomial 7(10) : 2 x0x7 +2 x1x6 +2 x2x10 +2 x2x11 +2 x3x10 -2 x3x11 +2 x4x8 +2 x4x9 +2 x5x8 -2 x5x9

	Polynomial 8(8) : 2 x0x8 -2 x1x9 +2 x2x6 -2 x3x7 +2 x4x10 -2 x4x11 -2 x5x10 -2 x5x11

	Polynomial 9(8) : 2 x0x9 +2 x1x8 +2 x2x7 +2 x3x6 +2 x4x10 +2 x4x11 +2 x5x10 -2 x5x11

	Polynomial 10(6) : 2 x0x10 -2 x1x11 +2 x2x8 -2 x3x9 +2 x4x6 -2 x5x7

	Polynomial 11(6) : 2 x0x11 +2 x1x10 +2 x2x9 +2 x3x8 +2 x4x7 +2 x5x6
	*/
	BLST_FORCEINLINE StandardElement fp12_flat_sqr(const StandardElement &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);
		// 0 1
		auto poly_0_pos_shared = Constants::offset() + ((x[2] * x[4] + x[3] * minus_x5 + x[6] * x[10] + x[7] * minus_x11) * _2) + x[8] * x[8] + x[9] * minus_x9;
		auto poly_0_neg_shared = (x[2] * x[5] + x[3] * x[4] + x[6] * x[11] + x[7] * x[10] + x[8] * x[9]) * _2;
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * _2;
		// 2 3
		auto poly_2_pos_shared = Constants::offset() + ((x[8] * x[10] + x[9] * minus_x11) * _2) + x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_neg_shared = (x[4] * x[5] + x[8] * x[11] + x[9] * x[10]) * _2;
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + ((x[0] * x[2] + x[1] * minus_x3) * _2) + x[6] * x[6] + x[7] * minus_x7;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + (x[0] * x[3] + x[1] * x[2] + x[6] * x[7]) * _2;
		// 4 5
		auto poly_4_pos_shared = x[10] * x[10] + x[11] * minus_x11;
		auto poly_4_unique = poly_4_pos_shared + ((x[0] * x[4] + x[1] * minus_x5 + x[6] * x[8] + x[7] * minus_x9 + x[10] * minus_x11) * _2) + x[2] * x[2] + x[3] * minus_x3;
		auto poly_5_unique = poly_4_pos_shared + (x[0] * x[5] + x[1] * x[4] + x[2] * x[3] + x[6] * x[9] + x[7] * x[8] + x[10] * x[11]) * _2;
		// 6 7
		auto poly_6_pos_shared = Constants::offset() + x[2] * x[10] + x[3] * minus_x11 + x[4] * x[8] + x[5] * minus_x9;
		auto poly_6_neg_shared = x[2] * x[11] + x[3] * x[10] + x[4] * x[9] + x[5] * x[8];
		auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[0] * x[6] + x[1] * minus_x7) * _2;
		auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[0] * x[7] + x[1] * x[6]) * _2;
		// 8 9
		auto poly_8_pos_shared = Constants::offset() + x[4] * x[10] + x[5] * minus_x11;
		auto poly_8_neg_shared = x[4] * x[11] + x[5] * x[10];
		auto poly_8_unique = (poly_8_pos_shared - poly_8_neg_shared + x[0] * x[8] + x[1] * minus_x9 + x[2] * x[6] + x[3] * minus_x7) * _2;
		auto poly_9_unique = (poly_8_pos_shared + poly_8_neg_shared + x[0] * x[9] + x[1] * x[8] + x[2] * x[7] + x[3] * x[6]) * _2;
		// 10 11
		auto poly_10_unique = (x[0] * x[10] + x[1] * minus_x11 + x[2] * x[8] + x[3] * minus_x9 + x[4] * x[6] + x[5] * minus_x7) * _2;
		auto poly_11_unique = (x[0] * x[11] + x[1] * x[10] + x[2] * x[9] + x[3] * x[8] + x[4] * x[7] + x[5] * x[6]) * _2;
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
	}
	/* 
	Polynomial 0(6) : 3 x0x0 -2 x0 -3 x1x1 +3 x8x8 -3 x9x9 -6 x8x9

	Polynomial 1(5) : -2 x1 +3 x8x8 -3 x9x9 +6 x0x1 +6 x8x9

	Polynomial 2(6) : -2 x2 +3 x4x4 -3 x5x5 +3 x6x6 -3 x7x7 -6 x4x5

	Polynomial 3(5) : -2 x3 +3 x4x4 -3 x5x5 +6 x4x5 +6 x6x7

	Polynomial 4(6) : 3 x2x2 -3 x3x3 -2 x4 +3 x10x10 -3 x11x11 -6 x10x11

	Polynomial 5(5) : -2 x5 +3 x10x10 -3 x11x11 +6 x2x3 +6 x10x11

	Polynomial 6(5) : 2 x6 +6 x2x10 -6 x2x11 -6 x3x10 -6 x3x11

	Polynomial 7(5) : 2 x7 +6 x2x10 +6 x2x11 +6 x3x10 -6 x3x11

	Polynomial 8(3) : 2 x8 +6 x0x8 -6 x1x9

	Polynomial 9(3) : 2 x9 +6 x0x9 +6 x1x8

	Polynomial 10(3) : 2 x10 +6 x4x6 -6 x5x7

	Polynomial 11(3) : 2 x11 +6 x4x7 +6 x5x6
	*/

	template<int64_t N>
	inline auto fp12_flat_cyclotomic_sqr(const CommonElement<N> &x, const Ring &ring) const {
		// Cyclotomic membership is a precondition. This fixed split is selected
		// because the monolithic form exhausts the register file. It performs no
		// validation, allocation, exception handling, or data-dependent dispatch.
		// Both noinline halves write into this final return owner. Keeping the input
		// distinct from that owner preserves exact caller-side storage reuse without
		// making the second half observe coordinates written by the first half.
		CommonElement<6> result;
		fp12_flat_cyclotomic_sqr_first_split(result, x, ring);
		fp12_flat_cyclotomic_sqr_second_split(result, x, ring);
		return result;
	}

	template<int64_t N>
	BLST_AVX_FP12_NOINLINE void fp12_flat_cyclotomic_sqr_first_split(
		CommonElement<6> &result, const CommonElement<N> &x, const Ring &ring) const {
		// Keep this source text aligned with the first accepted split proof ledger.
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);

		auto poly_0_pos_shared = x[8] * x[8] + x[9] * minus_x9;
		auto poly_0_unique = (poly_0_pos_shared + ((x[8] * minus_x9) * _2) + x[0] * x[0] + x[0] * Constants::neg_two_thirds() + x[1] * minus_x1);
		auto poly_1_unique = (poly_0_pos_shared + ((x[0] * x[1] + x[8] * x[9]) * _2) + x[1] * Constants::neg_two_thirds());
		auto poly_2_pos_shared = x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_unique = (poly_2_pos_shared + ((x[4] * minus_x5) * _2) + x[2] * Constants::neg_two_thirds() + x[6] * x[6] + x[7] * minus_x7);
		auto poly_3_unique = (poly_2_pos_shared + ((x[6] * x[7] + x[4] * x[5]) * _2) + x[3] * Constants::neg_two_thirds());
		auto poly_4_pos_shared = x[10] * x[10] + x[11] * minus_x11;
		auto poly_4_unique = (poly_4_pos_shared + ((x[10] * minus_x11) * _2) + x[2] * x[2] + x[3] * minus_x3 + x[4] * Constants::neg_two_thirds());
		auto poly_5_unique = (poly_4_pos_shared + ((x[2] * x[3] + x[10] * x[11]) * _2) + x[5] * Constants::neg_two_thirds());

		auto reduced = batch_reduce<first_cyclotomic_reduce_bound>(
			ring, poly_0_unique, poly_1_unique,
			poly_2_unique, poly_3_unique, poly_4_unique, poly_5_unique);
		auto scaled = to_common<6>(ring, reduced[0] * _3,
			reduced[1] * _3, reduced[2] * _3, reduced[3] * _3,
			reduced[4] * _3, reduced[5] * _3);
		const auto expanded = ring.template batch_expand<6>(scaled);
		result[0] = expanded[0];
		result[1] = expanded[1];
		result[2] = expanded[2];
		result[3] = expanded[3];
		result[4] = expanded[4];
		result[5] = expanded[5];
	}

	template<int64_t N>
	BLST_AVX_FP12_NOINLINE void fp12_flat_cyclotomic_sqr_second_split(
		CommonElement<6> &result, const CommonElement<N> &x, const Ring &ring) const {
		// Keep this source text aligned with the second accepted split proof ledger.
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);

		auto poly_6_pos_shared = (Constants::offset() + x[2] * x[10] + x[3] * minus_x11);
		auto poly_6_neg_shared = (x[2] * x[11] + x[3] * x[10]);
		auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[6] * Constants::one_third());
		auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[7] * Constants::one_third());
		auto poly_8_unique = (((x[0] * x[8] + x[1] * minus_x9)) + x[8] * Constants::one_third());
		auto poly_9_unique = (((x[0] * x[9] + x[1] * x[8])) + x[9] * Constants::one_third());
		auto poly_10_unique = (((x[4] * x[6] + x[5] * minus_x7)) + x[10] * Constants::one_third());
		auto poly_11_unique = (((x[4] * x[7] + x[5] * x[6])) + x[11] * Constants::one_third());

		auto reduced = batch_reduce<second_cyclotomic_reduce_bound>(
			ring, poly_6_unique, poly_7_unique,
			poly_8_unique, poly_9_unique, poly_10_unique, poly_11_unique);
		auto scaled = to_common<6>(ring, reduced[0] * _6,
			reduced[1] * _6, reduced[2] * _6, reduced[3] * _6,
			reduced[4] * _6, reduced[5] * _6);
		const auto expanded = ring.template batch_expand<6>(scaled);
		result[6] = expanded[0];
		result[7] = expanded[1];
		result[8] = expanded[2];
		result[9] = expanded[3];
		result[10] = expanded[4];
		result[11] = expanded[5];
	}

	inline StandardElement fp12_flat_frobenius_map1(const StandardElement &x, const Ring &ring) const {
		auto frob1_2_2 = x[2] * Constants::frob1_2();
		auto frob1_2_3 = x[3] * Constants::frob1_2();
		
		auto frob1_3_4 = x[4] * Constants::frob1_3();
		
		auto frob1_4_5 = x[5] * Constants::frob1_4();
		
		auto frob1_5_6 = x[6] * Constants::frob1_5();
		
		auto frob1_6_6 = x[6] * Constants::frob1_6();
		auto frob1_6_7 = x[7] * Constants::frob1_6();
		
		auto frob1_7_8 = x[8] * Constants::frob1_7();
		auto frob1_7_9 = x[9] * Constants::frob1_7();
		
		auto frob1_8_9 = x[9] * Constants::frob1_8();
		
		auto frob1_9_10 = x[10] * Constants::frob1_9();
		
		auto frob1_10_10 = x[10] * Constants::frob1_10();
		auto frob1_10_11 = x[11] * Constants::frob1_10();
		auto poly2_unique = frob1_2_3;
		auto poly3_unique = frob1_2_2;
		auto poly4_unique = frob1_3_4;
		auto poly5_unique = frob1_4_5;
		auto shared_frob1_6_7 = frob1_6_7;
		auto shared_frob1_7_8 = frob1_7_8;
		auto shared_frob1_10_11 = frob1_10_11;
		auto poly6_unique = shared_frob1_6_7 + frob1_5_6;
		auto poly7_unique = shared_frob1_6_7 + frob1_6_6;
		auto poly8_unique = shared_frob1_7_8 + frob1_7_9;
		auto poly9_unique = shared_frob1_7_8 + frob1_8_9;
		auto poly10_unique = shared_frob1_10_11 + frob1_9_10;
		auto poly11_unique = shared_frob1_10_11 + frob1_10_10;
		auto [reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11] = ring.batch_reduce_expand(poly2_unique, poly3_unique, poly4_unique, poly5_unique, poly6_unique, poly7_unique, poly8_unique, poly9_unique, poly10_unique, poly11_unique);
		return StandardElement{x[0], ring.standard_negate(x[1]), reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11};
	}
	
	template<int64_t N>
	inline auto fp12_flat_frobenius_map2(const CommonElement<N> &x, const Ring &ring) const {
		auto frob2_1_2 = x[2] * Constants::frob2_1();
		auto frob2_1_3 = x[3] * Constants::frob2_1();
		
		auto frob2_2_4 = x[4] * Constants::frob2_2();
		auto frob2_2_5 = x[5] * Constants::frob2_2();
		
		auto frob2_3_6 = x[6] * Constants::frob2_3();
		auto frob2_3_7 = x[7] * Constants::frob2_3();
		
		auto frob2_5_10 = x[10] * Constants::frob2_5();
		auto frob2_5_11 = x[11] * Constants::frob2_5();
		auto poly2_unique = frob2_1_2;
		auto poly3_unique = frob2_1_3;
		auto poly4_unique = frob2_2_4;
		auto poly5_unique = frob2_2_5;
		auto poly6_unique = frob2_3_6;
		auto poly7_unique = frob2_3_7;
		auto poly10_unique = frob2_5_10;
		auto poly11_unique = frob2_5_11;
		auto [reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, reduced_poly10, reduced_poly11] = ring.batch_reduce_expand(poly2_unique, poly3_unique, poly4_unique, poly5_unique, poly6_unique, poly7_unique, poly10_unique, poly11_unique);
		return CommonElement<N> {x[0].template recast_rns_bounds<0, N>(), x[1].template recast_rns_bounds<0, N>(), reduced_poly2.template recast_rns_bounds<0, N>(), reduced_poly3.template recast_rns_bounds<0, N>(), reduced_poly4.template recast_rns_bounds<0, N>(), reduced_poly5.template recast_rns_bounds<0, N>(), reduced_poly6.template recast_rns_bounds<0, N>(), reduced_poly7.template recast_rns_bounds<0, N>(), ring.standard_negate(x[8]).template recast_rns_bounds<0, N>(), ring.standard_negate(x[9]).template recast_rns_bounds<0, N>(), reduced_poly10.template recast_rns_bounds<0, N>(), reduced_poly11.template recast_rns_bounds<0, N>()};
	}
	
	inline StandardElement fp12_flat_frobenius_map3(const StandardElement &x, const Ring &ring) const {
		auto frob3_2_6 = x[6] * Constants::frob3_2();
		auto frob3_2_8 = x[8] * Constants::frob3_2();
		auto frob3_2_9 = x[9] * Constants::frob3_2();
		auto frob3_2_10 = x[10] * Constants::frob3_2();
		auto frob3_2_11 = x[11] * Constants::frob3_2();
		
		auto frob3_3_9 = x[9] * Constants::frob3_3();
		auto frob3_3_10 = x[10] * Constants::frob3_3();
		auto frob3_3_6 = x[6] * Constants::frob3_3();
		auto frob3_3_7 = x[7] * Constants::frob3_3();
		auto shared_frob3_3_7 = frob3_3_7;
		auto shared_frob3_2_8 = frob3_2_8;
		auto shared_frob3_2_11 = frob3_2_11;
		auto poly6_unique = shared_frob3_3_7 + frob3_2_6;
		auto poly7_unique = shared_frob3_3_7 + frob3_3_6;
		auto poly8_unique = shared_frob3_2_8 + frob3_2_9;
		auto poly9_unique = shared_frob3_2_8 + frob3_3_9;
		auto poly10_unique = shared_frob3_2_11 + frob3_3_10;
		auto poly11_unique = shared_frob3_2_11 + frob3_2_10;
		auto [reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11] = ring.batch_reduce_expand(poly6_unique, poly7_unique, poly8_unique, poly9_unique, poly10_unique, poly11_unique);
		return StandardElement{x[0], ring.standard_negate(x[1]), x[3], x[2], ring.standard_negate(x[4]), x[5], reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11};
	}
	
	template<int64_t N>
	inline auto fp12_flat_conjugate(const CommonElement<N> &x, const Ring &ring) const {
		return CommonElement<N> {
			x[0], x[1], x[2], x[3], x[4], x[5], ring.standard_negate(x[6]), ring.standard_negate(x[7]), ring.standard_negate(x[8]), ring.standard_negate(x[9]), ring.standard_negate(x[10]), ring.standard_negate(x[11])
		};
	}

	/* 
	Polynomial 0(13) : 1 x0x0 -1 x1x1 -1 x8x8 +1 x9x9 +2 x2x4 -2 x2x5 -2 x3x4 -2 x3x5 -2 x6x10 +2 x6x11 +2 x7x10 +2 x7x11 +2 x8x9

	Polynomial 1(12) : -1 x8x8 +1 x9x9 +2 x0x1 +2 x2x4 +2 x2x5 +2 x3x4 -2 x3x5 -2 x6x10 -2 x6x11 -2 x7x10 +2 x7x11 -2 x8x9

	Polynomial 2(11) : 1 x4x4 -1 x5x5 -1 x6x6 +1 x7x7 +2 x0x2 -2 x1x3 -2 x4x5 -2 x8x10 +2 x8x11 +2 x9x10 +2 x9x11

	Polynomial 3(10) : 1 x4x4 -1 x5x5 +2 x0x3 +2 x1x2 +2 x4x5 -2 x6x7 -2 x8x10 -2 x8x11 -2 x9x10 +2 x9x11

	Polynomial 4(9) : 1 x2x2 -1 x3x3 -1 x10x10 +1 x11x11 +2 x0x4 -2 x1x5 -2 x6x8 +2 x7x9 +2 x10x11

	Polynomial 5(8) : -1 x10x10 +1 x11x11 +2 x0x5 +2 x1x4 +2 x2x3 -2 x6x9 -2 x7x8 -2 x10x11
	*/
	inline auto fp12_flat_inverse_to_fp6(const std::array<typename Ring::StandardElement, 12> &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x6 = ring.negate(x[6]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x8 = ring.negate(x[8]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x10 = ring.negate(x[10]);
		auto minus_x11 = ring.negate(x[11]);
	// 0 1
	auto poly_0_pos_shared = Constants::offset() + ((x[2] * x[4] + x[3] * minus_x5 + x[6] * minus_x10 + x[7] * x[11]) * _2) + x[8] * minus_x8 + x[9] * x[9];
	auto poly_0_neg_shared = (x[2] * x[5] + x[3] * x[4] + x[6] * minus_x11 + x[7] * minus_x10 + x[8] * minus_x9) * _2;
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * _2;
	// 2 3
	auto poly_2_pos_shared = Constants::offset() + ((x[8] * minus_x10 + x[9] * x[11]) * _2) + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_neg_shared = (x[4] * x[5] + x[8] * minus_x11 + x[9] * minus_x10) * _2;
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + ((x[0] * x[2] + x[1] * minus_x3) * _2) + x[6] * minus_x6 + x[7] * x[7];
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + (x[0] * x[3] + x[1] * x[2] + x[6] * minus_x7) * _2;
	// 4 5
	auto poly_4_pos_shared = Constants::offset() + x[10] * minus_x10 + x[11] * x[11];
	auto poly_4_unique = poly_4_pos_shared + ((x[0] * x[4] + x[1] * minus_x5 + x[6] * minus_x8 + x[7] * x[9] + x[10] * x[11]) * _2) + x[2] * x[2] + x[3] * minus_x3;
	auto poly_5_unique = poly_4_pos_shared + (x[0] * x[5] + x[1] * x[4] + x[2] * x[3] + x[6] * minus_x9 + x[7] * minus_x8 + x[10] * minus_x11) * _2;
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
	}
	/* 
	Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3

	Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3

	Polynomial 2(8) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5

	Polynomial 3(8) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5

	Polynomial 4(6) : 1 x0y4 -1 x1y5 +1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1

	Polynomial 5(6) : 1 x0y5 +1 x1y4 +1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0

	Polynomial 6(10) : -1 x6y0 +1 x7y1 -1 x8y4 +1 x8y5 +1 x9y4 +1 x9y5 -1 x10y2 +1 x10y3 +1 x11y2 +1 x11y3

	Polynomial 7(10) : -1 x6y1 -1 x7y0 -1 x8y4 -1 x8y5 -1 x9y4 +1 x9y5 -1 x10y2 -1 x10y3 -1 x11y2 +1 x11y3

	Polynomial 8(8) : -1 x6y2 +1 x7y3 -1 x8y0 +1 x9y1 -1 x10y4 +1 x10y5 +1 x11y4 +1 x11y5

	Polynomial 9(8) : -1 x6y3 -1 x7y2 -1 x8y1 -1 x9y0 -1 x10y4 -1 x10y5 -1 x11y4 +1 x11y5

	Polynomial 10(6) : -1 x6y4 +1 x7y5 -1 x8y2 +1 x9y3 -1 x10y0 +1 x11y1

	Polynomial 11(6) : -1 x6y5 -1 x7y4 -1 x8y3 -1 x9y2 -1 x10y1 -1 x11y0
	*/
	inline auto fp6_flat_inverse_to_fp12(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 6> &y, const Ring &ring) const {
		auto minus_y0 = ring.negate(y[0]);
		auto minus_y1 = ring.negate(y[1]);
		auto minus_y2 = ring.negate(y[2]);
		auto minus_y3 = ring.negate(y[3]);
		auto minus_y4 = ring.negate(y[4]);
		auto minus_y5 = ring.negate(y[5]);
		// 0 1
		auto poly_0_pos_shared = Constants::offset() + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3;
		auto poly_0_neg_shared = x[2] * y[5] + x[3] * y[4] + x[4] * y[3] + x[5] * y[2];
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
		// 2 3
		auto poly_2_pos_shared = Constants::offset() + x[4] * y[4] + x[5] * minus_y5;
		auto poly_2_neg_shared = x[4] * y[5] + x[5] * y[4];
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y[2] + x[1] * minus_y3 + x[2] * y[0] + x[3] * minus_y1;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0];
		// 4 5
		auto poly_4_unique = x[0] * y[4] + x[1] * minus_y5 + x[2] * y[2] + x[3] * minus_y3 + x[4] * y[0] + x[5] * minus_y1;
		auto poly_5_unique = x[0] * y[5] + x[1] * y[4] + x[2] * y[3] + x[3] * y[2] + x[4] * y[1] + x[5] * y[0];
		// 6 7
		auto poly_6_pos_shared = Constants::offset() + x[8] * minus_y4 + x[9] * y[5] + x[10] * minus_y2 + x[11] * y[3];
		auto poly_6_neg_shared = x[8] * minus_y5 + x[9] * minus_y4 + x[10] * minus_y3 + x[11] * minus_y2;
		auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[6] * minus_y0 + x[7] * y[1];
		auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[6] * minus_y1 + x[7] * minus_y0;
		// 8 9
		auto poly_8_pos_shared = Constants::offset() + x[10] * minus_y4 + x[11] * y[5];
		auto poly_8_neg_shared = x[10] * minus_y5 + x[11] * minus_y4;
		auto poly_8_unique = poly_8_pos_shared - poly_8_neg_shared + x[6] * minus_y2 + x[7] * y[3] + x[8] * minus_y0 + x[9] * y[1];
		auto poly_9_unique = poly_8_pos_shared + poly_8_neg_shared + x[6] * minus_y3 + x[7] * minus_y2 + x[8] * minus_y1 + x[9] * minus_y0;
		// 10 11
		auto poly_10_unique = x[6] * minus_y4 + x[7] * y[5] + x[8] * minus_y2 + x[9] * y[3] + x[10] * minus_y0 + x[11] * y[1];
		auto poly_11_unique = x[6] * minus_y5 + x[7] * minus_y4 + x[8] * minus_y3 + x[9] * minus_y2 + x[10] * minus_y1 + x[11] * minus_y0;
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
	}
	/* 
	Polynomial 0(6) : 1 x0x0 -1 x1x1 -1 x2x4 +1 x2x5 +1 x3x4 +1 x3x5

	Polynomial 1(5) : 2 x0x1 -1 x2x4 -1 x2x5 -1 x3x4 +1 x3x5

	Polynomial 2(5) : 1 x4x4 -1 x5x5 -1 x0x2 +1 x1x3 -2 x4x5

	Polynomial 3(5) : 1 x4x4 -1 x5x5 -1 x0x3 -1 x1x2 +2 x4x5

	Polynomial 4(4) : 1 x2x2 -1 x3x3 -1 x0x4 +1 x1x5

	Polynomial 5(3) : -1 x0x5 -1 x1x4 +2 x2x3
	*/
	inline auto fp6_flat_inverse_to_fp2_c0_c1_c2(const std::array<typename Ring::StandardElement, 6> &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x2 = ring.negate(x[2]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x4 = ring.negate(x[4]);
		auto minus_x5 = ring.negate(x[5]);
		// 0 1
		auto poly_0_pos_shared = Constants::offset() + x[2] * minus_x4 + x[3] * x[5];
		auto poly_0_neg_shared = x[2] * minus_x5 + x[3] * minus_x4;
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * _2;
		// 2 3
		auto poly_2_pos_shared = Constants::offset() + x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_unique = poly_2_pos_shared + ((x[4] * minus_x5) * _2) + x[0] * minus_x2 + x[1] * x[3];
		auto poly_3_unique = poly_2_pos_shared + ((x[4] * x[5]) * _2) + x[0] * minus_x3 + x[1] * minus_x2;
		// 4 5
		auto poly_4_unique = x[2] * x[2] + x[3] * minus_x3 + x[0] * minus_x4 + x[1] * x[5];
		auto poly_5_unique = ((x[2] * x[3]) * _2) + x[0] * minus_x5 + x[1] * minus_x4;
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
	}
	/* 
	Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3

	Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3
	*/
	inline auto fp6_flat_inverse_to_fp2_ret(const std::array<typename Ring::StandardElement, 6> &x, const std::array<typename Ring::StandardElement, 6> &y, const Ring &ring) const {
		auto minus_y1 = ring.negate(y[1]);
		auto minus_y3 = ring.negate(y[3]);
		auto minus_y5 = ring.negate(y[5]);
		// 0 1
		auto poly_0_pos_shared = Constants::offset() + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3;
		auto poly_0_neg_shared = x[2] * y[5] + x[3] * y[4] + x[4] * y[3] + x[5] * y[2];
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique);
	}
	/* 
	Polynomial 0(2) : 1 x0x2 -1 x1x3

	Polynomial 1(2) : 1 x0x3 +1 x1x2

	Polynomial 2(2) : 1 x0x4 -1 x1x5

	Polynomial 3(2) : 1 x0x5 +1 x1x4

	Polynomial 4(2) : 1 x0y0 -1 x1y1

	Polynomial 5(2) : 1 x0y1 +1 x1y0
	*/
	inline auto fp2_flat_inverse_to_fp6(const std::array<typename Ring::StandardElement, 6> &x, const std::array<typename Ring::StandardElement, 2> &y, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		// 0 1
		auto poly_0_unique = y[0] * x[0] + y[1] * minus_x1;
		auto poly_1_unique = y[0] * x[1] + y[1] * x[0];
		// 2 3
		auto poly_2_unique = y[0] * x[2] + y[1] * minus_x3;
		auto poly_3_unique = y[0] * x[3] + y[1] * x[2];
		// 4 5
		auto poly_4_unique = y[0] * x[4] + y[1] * minus_x5;
		auto poly_5_unique = y[0] * x[5] + y[1] * x[4];
		
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
	}
	/* 
	Polynomial 0(2) : 1 x0x0 +1 x1x1
	*/
	inline auto fp2_flat_inverse_to_fp(const std::array<typename Ring::StandardElement, 2> &x, const Ring &ring) const {
		auto poly_0_unique = x[0] * x[0] + x[1] * x[1];
		return ring.batch_reduce_expand(poly_0_unique);
	}
	/* 
	Polynomial 0(1) : 1 x0y0

	Polynomial 1(1) : -1 x1y0
	*/
	inline auto fp_flat_inverse_to_fp2(const std::array<typename Ring::StandardElement, 2> &x, const std::array<typename Ring::StandardElement, 1> &y, const Ring &ring) const {
		auto minus_y0 = ring.negate(y[0]);
		// 0 1
		auto poly_0_unique = x[0] * y[0];
		auto poly_1_unique = x[1] * minus_y0;
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique);
	}


	template<class Inverter>
	inline auto fp12_flat_inverse(const StandardElement &x, const Inverter &inverter, const Ring &ring) const {
		// Precondition: x != 0. No check is performed here. Timing and failure
		// behavior of the base-field reciprocal are exactly those of `inverter`.
		auto f6 = fp12_flat_inverse_to_fp6(x, ring);
		auto f2_c0_c1_c2 = fp6_flat_inverse_to_fp2_c0_c1_c2(f6, ring);
		auto f2 = fp6_flat_inverse_to_fp2_ret(f2_c0_c1_c2, f6, ring);
		auto f1 = fp2_flat_inverse_to_fp(f2, ring);
		const auto& f1_elem = f1[0];
		auto f1_inv = inverter.invert(f1_elem, ring);
		std::array<typename Ring::StandardElement, 1> f1_inv_arr = {f1_inv};
		auto f2_inv = fp_flat_inverse_to_fp2(f2, f1_inv_arr, ring);
		auto f6_inv = fp2_flat_inverse_to_fp6(f2_c0_c1_c2, f2_inv, ring);
		return fp6_flat_inverse_to_fp12(x, f6_inv, ring);
	}
};

}  // namespace blst::avx

#undef BLST_AVX_FP12_NOINLINE
