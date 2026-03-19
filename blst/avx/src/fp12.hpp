#pragma once
#include <array>
#include <cstdint>
#include <type_traits>
#include <iostream>
#include "../cpu/precompute/gmp_wrapper.hpp"
#include "fp2.hpp"

/*
template<class... Args>
class FP12Element {
	public:
	constexpr static int size = sizeof...(Args);
	std::tuple<Args...> elements;

	INLINE FP12Element(const Args &... args)
		: elements(args...) {}

	INLINE FP12Element()
		: elements(Args()...) {}

	template<int index>
	INLINE auto operator[](std::integral_constant<int, index>) const {
		return std::get<index>(elements);
	}

	template<int index>
	INLINE auto operator[](std::integral_constant<int, index>) {
		return std::get<index>(elements);
	}
};
*/

// concat is defined in multiplication.hpp

template<class Ring, int VERSION=3>
class FP12Ring {

	public:
	// using StandardElement = std::array<typename Ring::StandardElement, 12>;
	using FP2StandardElement = typename FP2Ring<Ring>::StandardElement;

	// For simplicity we recast up to the maximum and use an array
	// You could also use a tuple for more fine grained control.
	template<int64_t common_bound>
	using CommonElement = std::array<typename Ring::template ReducedElement<common_bound>, 12>;

	using StandardElement = CommonElement<1>;
	
	// Shorthand for integral constants
	static constexpr auto _2 = std::integral_constant<int, 2>{};
	static constexpr auto _3 = std::integral_constant<int, 3>{};
	static constexpr auto _6 = std::integral_constant<int, 6>{};

	private:
	typename Ring::StandardElement one_third;
	typename Ring::StandardElement neg_one_third;
	typename Ring::StandardElement two_thirds;
	typename Ring::StandardElement neg_two_thirds;
	typename Ring::StandardElement frob1_2;
	typename Ring::StandardElement frob1_3;
	typename Ring::StandardElement frob1_4;
	typename Ring::StandardElement frob1_5;
	typename Ring::StandardElement frob1_6;
	typename Ring::StandardElement frob1_7;
	typename Ring::StandardElement frob1_8;
	typename Ring::StandardElement frob1_9;
	typename Ring::StandardElement frob1_10;
	typename Ring::StandardElement frob2_1;
	typename Ring::StandardElement frob2_2;
	typename Ring::StandardElement frob2_3;
	typename Ring::StandardElement frob2_5;
	typename Ring::StandardElement frob3_2;
	typename Ring::StandardElement frob3_3;
	typename Ring::StandardElement frob3_4;

	// We do ten subtractions, so this is 10*RNS_modulus + 10*modulus to ensure everything remains positive.
	constexpr static int offset_rns_bound = 10 * Ring::STANDARD_BOUND * Ring::STANDARD_BOUND * (VERSION == 1 ? 1 : 6 * 6);
	constexpr static int offset_ele_bound = 10 * Ring::STANDARD_BOUND * std::max(Ring::STANDARD_BOUND, Ring::MULT_OK_BOUND);
	typename Ring::template WideOffsetElement<offset_ele_bound, offset_rns_bound> offset;

	public:
	FP12Ring(const Ring &ring) {
		// Initialize all of the constants needed for the flattened code
		one_third = ring.from_hex("0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71d");
		neg_one_third = ring.from_hex("0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38e");
		two_thirds = ring.from_hex("0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38f");
		neg_two_thirds = ring.from_hex("0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c");
		frob1_2 = ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac");
		frob1_3 = ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad");
		frob1_4 = ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe");
		frob1_5 = ring.from_hex("0x1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8");
		frob1_6 = ring.from_hex("0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3");
		frob1_7 = ring.from_hex("0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09");
		frob1_8 = ring.from_hex("0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2");
		frob1_9 = ring.from_hex("0x5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116");
		frob1_10 = ring.from_hex("0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995");
		frob2_1 = ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe");
		frob2_2 = frob1_2;
		frob2_3 = ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff");
		frob2_5 = ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad");
		frob3_2 = ring.from_hex("0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2");
		frob3_3 = ring.from_hex("0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09");
		offset = ring.template compute_wide_offset<offset_ele_bound, offset_rns_bound>();
	}
    /*
	auto expand(const CommonElement &x, const Ring &ring) const {
		return ring.template batch_expand<12>(x);
	}

	// no-op for already expanded elements
	template<class ElementType, class RNSBounds>
	auto expand(const std::array<ExpandedElement<ElementType, RNSBounds>, 12> &x, const Ring &ring) const {
		return x;
	}
	*/

	
	template<int64_t upper_bound, class E1, class E2, class E3, class E4, class E5, class E6, class E7, class E8, class E9, class E10, class E11, class E12>
	INLINE auto to_common(const Ring &ring, const E1 &e1, const E2 &e2, const E3 &e3, const E4 &e4, const E5 &e5, const E6 &e6, const E7 &e7, const E8 &e8, const E9 &e9, const E10 &e10, const E11 &e11, const E12 &e12) const {
		using ResultType = typename Ring::template SmallReducedElement<upper_bound>;
		return std::array<ResultType, 12>{{
			ring.prep_expand(e1).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e2).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e3).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e4).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e5).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e6).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e7).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e8).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e9).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e10).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e11).template recast_rns_bounds<0, upper_bound>(),
			ring.prep_expand(e12).template recast_rns_bounds<0, upper_bound>(),
		}};
	}

	template<int64_t upper_bound, class E1, class E2, class E3, class E4, class E5, class E6>
	INLINE auto to_common(const Ring &ring, const E1 &e1, const E2 &e2, const E3 &e3, const E4 &e4, const E5 &e5, const E6 &e6) const {
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
	INLINE auto recast_common(const CommonElement<Z> &x) const {
		constexpr int N = VERSION == 1 ? 1 : 6;
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

	template<typename... Args>
	INLINE auto batch_reduce(const Ring &ring, Args&&... args) const {
		constexpr int batch = sizeof...(Args);
		auto tuple = std::tie(args...);
		auto ready_array = make_ready_array_impl_22(tuple, std::make_index_sequence<batch>{}, ring);
		return ring.template batch_reduce<batch>(ready_array);
	}

	template<typename... Args, std::size_t... I>
	INLINE auto make_ready_array_impl_22(const std::tuple<Args...> &args, std::index_sequence<I...>, const Ring &ring) const {
		return std::array<typename Ring::template ReadyElement<200>, sizeof...(Args)>{ring.template ready<200>(std::get<I>(args))...};
	}

	StandardElement one(const Ring &ring) const {
		return std::array<typename Ring::StandardElement, 12> {ring.one(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), ring.zero()};
	}

	bool check_bounds(const Ring & ring) const {
		auto r = ring.check_bounds(offset, "fp12 offset");
		if (r.first != BigInt(0)) {
			std::cout << "Value was " << r.first.to_string(16) << " but expected 0" << std::endl;
		}
		return r.second;
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
	
/*
	template<bool recast = true, class Y>
	inline auto fp12_flat_mul(const CommonElement &x, const Y &y, const Ring &ring) const {
		// expand (acc = offset m1) if needed, do math (m1) on output, elementwise reduce, store
		auto m1_result = fp12_flat_mul_i<1>(x, y, ring);
		// reduce (acc = offset m2), do math (m2) on output, elementwise reduce, scalar multiply, store
		StandardElement m2_result = fp12_flat_mul_i<2>(x, y, ring, m1_result);
		if constexpr (recast) {
			return CommonElement(m2_result);
		} else {
			return m2_result;
		}
	}

	template<int mod_idx, class Y>
	inline StandardElement fp12_flat_mul_i(const CommonElement &x, const Y &y_in, const Ring &ring, const StandardElement&m1_result=StandardElement{}, const RNSMatrix2 &mat, const MontgomeryReduce2 &reducer) const {
		WideCommonElement y;
		if constexpr (mod_idx == 2) || (!expanded) {
			y = matmul(y_in, ring, mat);
		} else {
			y = M1(y_in);
		}
		auto minus_y1 = ring.negate(y[1], reducer);
		auto minus_y3 = ring.negate(y[3], reducer);
		auto minus_y5 = ring.negate(y[5], reducer);
		auto minus_y7 = ring.negate(y[7], reducer);
		auto minus_y9 = ring.negate(y[9], reducer);
		auto minus_y11 = ring.negate(y[11], reducer);

		// Then remove offset and do +=?
*/

	template<int64_t N1, int64_t N2>
	inline auto fp12_flat_mul(const CommonElement<N1> &x, const CommonElement<N2> &y, const Ring &ring) const {
		auto minus_y1 = ring.negate(y[1]);
		auto minus_y3 = ring.negate(y[3]);
		auto minus_y5 = ring.negate(y[5]);
		auto minus_y7 = ring.negate(y[7]);
		auto minus_y9 = ring.negate(y[9]);
		auto minus_y11 = ring.negate(y[11]);
		// 0 1
		auto poly_0_pos_shared = offset + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3 + x[6] * y[10] + x[7] * minus_y11 + x[8] * y[8] + x[9] * minus_y9 + x[10] * y[6] + x[11] * minus_y7;
		auto poly_0_neg_shared = x[2] * y[5] + x[3] * y[4] + x[4] * y[3] + x[5] * y[2] + x[6] * y[11] + x[7] * y[10] + x[8] * y[9] + x[9] * y[8] + x[10] * y[7] + x[11] * y[6];
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
		// 2 3
		auto poly_2_pos_shared = offset + x[4] * y[4] + x[5] * minus_y5 + x[8] * y[10] + x[9] * minus_y11 + x[10] * y[8] + x[11] * minus_y9;
		auto poly_2_neg_shared = x[4] * y[5] + x[5] * y[4] + x[8] * y[11] + x[9] * y[10] + x[10] * y[9] + x[11] * y[8];
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y[2] + x[1] * minus_y3 + x[2] * y[0] + x[3] * minus_y1 + x[6] * y[6] + x[7] * minus_y7;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0] + x[6] * y[7] + x[7] * y[6];
		// 4 5
		auto poly_4_pos_shared = offset + x[10] * y[10] + x[11] * minus_y11;
		auto poly_4_neg_shared = x[10] * y[11] + x[11] * y[10];
		auto poly_4_unique = poly_4_pos_shared - poly_4_neg_shared + x[0] * y[4] + x[1] * minus_y5 + x[2] * y[2] + x[3] * minus_y3 + x[4] * y[0] + x[5] * minus_y1 + x[6] * y[8] + x[7] * minus_y9 + x[8] * y[6] + x[9] * minus_y7;
		auto poly_5_unique = poly_4_pos_shared + poly_4_neg_shared + x[0] * y[5] + x[1] * y[4] + x[2] * y[3] + x[3] * y[2] + x[4] * y[1] + x[5] * y[0] + x[6] * y[9] + x[7] * y[8] + x[8] * y[7] + x[9] * y[6];
		// 6 7
		auto poly_6_pos_shared = offset + x[2] * y[10] + x[3] * minus_y11 + x[4] * y[8] + x[5] * minus_y9 + x[8] * y[4] + x[9] * minus_y5 + x[10] * y[2] + x[11] * minus_y3;
		auto poly_6_neg_shared = x[2] * y[11] + x[3] * y[10] + x[4] * y[9] + x[5] * y[8] + x[8] * y[5] + x[9] * y[4] + x[10] * y[3] + x[11] * y[2];
		auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[0] * y[6] + x[1] * minus_y7 + x[6] * y[0] + x[7] * minus_y1;
		auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[0] * y[7] + x[1] * y[6] + x[6] * y[1] + x[7] * y[0];
		// 8 9
		auto poly_8_pos_shared = offset + x[4] * y[10] + x[5] * minus_y11 + x[10] * y[4] + x[11] * minus_y5;
		auto poly_8_neg_shared = x[4] * y[11] + x[5] * y[10] + x[10] * y[5] + x[11] * y[4];
		auto poly_8_unique = poly_8_pos_shared - poly_8_neg_shared + x[0] * y[8] + x[1] * minus_y9 + x[2] * y[6] + x[3] * minus_y7 + x[6] * y[2] + x[7] * minus_y3 + x[8] * y[0] + x[9] * minus_y1;
		auto poly_9_unique = poly_8_pos_shared + poly_8_neg_shared + x[0] * y[9] + x[1] * y[8] + x[2] * y[7] + x[3] * y[6] + x[6] * y[3] + x[7] * y[2] + x[8] * y[1] + x[9] * y[0];
		// 10 11
		auto poly_10_unique = x[0] * y[10] + x[1] * minus_y11 + x[2] * y[8] + x[3] * minus_y9 + x[4] * y[6] + x[5] * minus_y7 + x[6] * y[4] + x[7] * minus_y5 + x[8] * y[2] + x[9] * minus_y3 + x[10] * y[0] + x[11] * minus_y1;
		auto poly_11_unique = x[0] * y[11] + x[1] * y[10] + x[2] * y[9] + x[3] * y[8] + x[4] * y[7] + x[5] * y[6] + x[6] * y[5] + x[7] * y[4] + x[8] * y[3] + x[9] * y[2] + x[10] * y[1] + x[11] * y[0];
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
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
        INLINE void flat_mulxy00z0(StandardElement &result, const StandardElement &x, const FP2<CX, CY> &constant_coord, const FP2StandardElement &x_coord, const FP2StandardElement &y_coord, const Ring &ring) const {
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
		auto poly_0_pos_shared = offset + x[4] * y2 + x[5] * minus_y3 + x[8] * y8 + x[9] * minus_y9;
		auto poly_0_neg_shared = x[4] * y3 + x[5] * y2 + x[8] * y9 + x[9] * y8;
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y0 + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y1 + x[1] * y0;
		// 2 3
		auto poly_2_pos_shared = offset + x[10] * y8 + x[11] * minus_y9;
		auto poly_2_neg_shared = x[10] * y9 + x[11] * y8;
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y2 + x[1] * minus_y3 + x[2] * y0 + x[3] * minus_y1;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y3 + x[1] * y2 + x[2] * y1 + x[3] * y0;
		// 4 5
		auto poly_4_unique = x[2] * y2 + x[3] * minus_y3 + x[4] * y0 + x[5] * minus_y1 + x[6] * y8 + x[7] * minus_y9;
		auto poly_5_unique = x[2] * y3 + x[3] * y2 + x[4] * y1 + x[5] * y0 + x[6] * y9 + x[7] * y8;
		// 6 7
		auto poly_6_pos_shared = offset + x[4] * y8 + x[5] * minus_y9 + x[10] * y2 + x[11] * minus_y3;
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
	INLINE StandardElement fp12_flat_sqr(const StandardElement &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);
		// 0 1
		auto poly_0_pos_shared = offset + ((x[2] * x[4] + x[3] * minus_x5 + x[6] * x[10] + x[7] * minus_x11) * _2) + x[8] * x[8] + x[9] * minus_x9;
		auto poly_0_neg_shared = (x[2] * x[5] + x[3] * x[4] + x[6] * x[11] + x[7] * x[10] + x[8] * x[9]) * _2;
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * _2;
		// 2 3
		auto poly_2_pos_shared = offset + ((x[8] * x[10] + x[9] * minus_x11) * _2) + x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_neg_shared = (x[4] * x[5] + x[8] * x[11] + x[9] * x[10]) * _2;
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + ((x[0] * x[2] + x[1] * minus_x3) * _2) + x[6] * x[6] + x[7] * minus_x7;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + (x[0] * x[3] + x[1] * x[2] + x[6] * x[7]) * _2;
		// 4 5
		auto poly_4_pos_shared = x[10] * x[10] + x[11] * minus_x11;
		auto poly_4_unique = poly_4_pos_shared + ((x[0] * x[4] + x[1] * minus_x5 + x[6] * x[8] + x[7] * minus_x9 + x[10] * minus_x11) * _2) + x[2] * x[2] + x[3] * minus_x3;
		auto poly_5_unique = poly_4_pos_shared + (x[0] * x[5] + x[1] * x[4] + x[2] * x[3] + x[6] * x[9] + x[7] * x[8] + x[10] * x[11]) * _2;
		// 6 7
		auto poly_6_pos_shared = offset + x[2] * x[10] + x[3] * minus_x11 + x[4] * x[8] + x[5] * minus_x9;
		auto poly_6_neg_shared = x[2] * x[11] + x[3] * x[10] + x[4] * x[9] + x[5] * x[8];
		auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[0] * x[6] + x[1] * minus_x7) * _2;
		auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[0] * x[7] + x[1] * x[6]) * _2;
		// 8 9
		auto poly_8_pos_shared = offset + x[4] * x[10] + x[5] * minus_x11;
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
		if constexpr (VERSION == 1) {
			// Old version, in case of regression
			return fp12_flat_cyclotomic_sqr_v1(x, ring);
		} else if constexpr (VERSION == 2) {
			// Batch of 12, multiply half by 3 and half by 6
			return fp12_flat_cyclotomic_sqr_v2(x, ring);
		} else {
			// Batch of 6 times 3, Batch of 6 times 6
			return fp12_flat_cyclotomic_sqr_v3(x, ring);
		}
	}

		
	template<int64_t N>
	inline auto fp12_flat_cyclotomic_sqr_v3(const CommonElement<N> &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);
		
		// 0 1
		auto poly_0_pos_shared = x[8] * x[8] + x[9] * minus_x9;
		auto poly_0_unique = (poly_0_pos_shared + ((x[8] * minus_x9) * _2) + x[0] * x[0] + x[0] * neg_two_thirds + x[1] * minus_x1);
		auto poly_1_unique = (poly_0_pos_shared + ((x[0] * x[1] + x[8] * x[9]) * _2) + x[1] * neg_two_thirds);
		// 2 3
		auto poly_2_pos_shared = x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_unique = (poly_2_pos_shared + ((x[4] * minus_x5) * _2) + x[2] * neg_two_thirds + x[6] * x[6] + x[7] * minus_x7);
		auto poly_3_unique = (poly_2_pos_shared + ((x[6] * x[7] + x[4] * x[5]) * _2) + x[3] * neg_two_thirds);
		// 4 5
		auto poly_4_pos_shared = x[10] * x[10] + x[11] * minus_x11;
		auto poly_4_unique = (poly_4_pos_shared + ((x[10] * minus_x11) * _2) + x[2] * x[2] + x[3] * minus_x3 + x[4] * neg_two_thirds);
		auto poly_5_unique = (poly_4_pos_shared + ((x[2] * x[3] + x[10] * x[11]) * _2) + x[5] * neg_two_thirds);

		auto reduced_first_half = batch_reduce(ring, poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
		auto reduced_mult_3_cast_6 = to_common<6>(ring, reduced_first_half[0] * _3, reduced_first_half[1] * _3, reduced_first_half[2] * _3, reduced_first_half[3] * _3, reduced_first_half[4] * _3, reduced_first_half[5] * _3);
		auto expanded_first_half = ring.template batch_expand<6>(reduced_mult_3_cast_6);
		// 6 7
		auto poly_6_pos_shared = (offset + x[2] * x[10] + x[3] * minus_x11);
		auto poly_6_neg_shared = (x[2] * x[11] + x[3] * x[10]);
		auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[6] * one_third);
		auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[7] * one_third);
		// 8 9
		auto poly_8_unique = (((x[0] * x[8] + x[1] * minus_x9)) + x[8] * one_third);
		auto poly_9_unique = (((x[0] * x[9] + x[1] * x[8])) + x[9] * one_third);
		// 10 11
		auto poly_10_unique = (((x[4] * x[6] + x[5] * minus_x7)) + x[10] * one_third);
		auto poly_11_unique = (((x[4] * x[7] + x[5] * x[6])) + x[11] * one_third);
		auto reduced_second_half = batch_reduce(ring, poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
		auto reduced_mult_6_cast_3 = to_common<6>(ring, reduced_second_half[0] * _6, reduced_second_half[1] * _6, reduced_second_half[2] * _6, reduced_second_half[3] * _6, reduced_second_half[4] * _6, reduced_second_half[5] * _6);
		auto expanded_second_half = ring.template batch_expand<6>(reduced_mult_6_cast_3);
		return concat(expanded_first_half, expanded_second_half);
	}
	
	template<int64_t N>
	inline auto fp12_flat_cyclotomic_sqr_v2(const CommonElement<N> &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);
		
		// 0 1
		auto poly_0_pos_shared = x[8] * x[8] + x[9] * minus_x9;
		auto poly_0_unique = (poly_0_pos_shared + ((x[8] * minus_x9) * _2) + x[0] * x[0] + x[0] * neg_two_thirds + x[1] * minus_x1);
		auto poly_1_unique = (poly_0_pos_shared + ((x[0] * x[1] + x[8] * x[9]) * _2) + x[1] * neg_two_thirds);
		// 2 3
		auto poly_2_pos_shared = x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_unique = (poly_2_pos_shared + ((x[4] * minus_x5) * _2) + x[2] * neg_two_thirds + x[6] * x[6] + x[7] * minus_x7);
		auto poly_3_unique = (poly_2_pos_shared + ((x[6] * x[7] + x[4] * x[5]) * _2) + x[3] * neg_two_thirds);
		// 4 5
		auto poly_4_pos_shared = x[10] * x[10] + x[11] * minus_x11;
		auto poly_4_unique = (poly_4_pos_shared + ((x[10] * minus_x11) * _2) + x[2] * x[2] + x[3] * minus_x3 + x[4] * neg_two_thirds);
		auto poly_5_unique = (poly_4_pos_shared + ((x[2] * x[3] + x[10] * x[11]) * _2) + x[5] * neg_two_thirds);
		// 6 7
		auto poly_6_pos_shared = (offset + x[2] * x[10] + x[3] * minus_x11);
		auto poly_6_neg_shared = (x[2] * x[11] + x[3] * x[10]);
		auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[6] * one_third);
		auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[7] * one_third);
		// 8 9
		auto poly_8_unique = (((x[0] * x[8] + x[1] * minus_x9)) + x[8] * one_third);
		auto poly_9_unique = (((x[0] * x[9] + x[1] * x[8])) + x[9] * one_third);
		// 10 11
		auto poly_10_unique = (((x[4] * x[6] + x[5] * minus_x7)) + x[10] * one_third);
		auto poly_11_unique = (((x[4] * x[7] + x[5] * x[6])) + x[11] * one_third);
		auto reduced = batch_reduce(ring, poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
		auto reduced_mult_3 = to_common<6>(ring, reduced[0] * _3, reduced[1] * _3, reduced[2] * _3, reduced[3] * _3, reduced[4] * _3, reduced[5] * _3, reduced[6] * _6, reduced[7] * _6, reduced[8] * _6, reduced[9] * _6, reduced[10] * _6, reduced[11] * _6);
		return ring.template batch_expand<12>(reduced_mult_3);
	}

	template<int64_t N>
	inline auto fp12_flat_cyclotomic_sqr_v1(const CommonElement<N> &x, const Ring &ring) const {
		auto minus_x1 = ring.negate(x[1]);
		auto minus_x3 = ring.negate(x[3]);
		auto minus_x5 = ring.negate(x[5]);
		auto minus_x7 = ring.negate(x[7]);
		auto minus_x9 = ring.negate(x[9]);
		auto minus_x11 = ring.negate(x[11]);
		
		// 0 1
		auto poly_0_pos_shared = x[8] * x[8] + x[9] * minus_x9;
		auto poly_0_unique = (poly_0_pos_shared + ((x[8] * minus_x9) * _2) + x[0] * x[0] + x[0] * neg_two_thirds + x[1] * minus_x1) * _3;
		auto poly_1_unique = (poly_0_pos_shared + ((x[0] * x[1] + x[8] * x[9]) * _2) + x[1] * neg_two_thirds) * _3;
		// 2 3
		auto poly_2_pos_shared = x[4] * x[4] + x[5] * minus_x5;
		auto poly_2_unique = (poly_2_pos_shared + ((x[4] * minus_x5) * _2) + x[2] * neg_two_thirds + x[6] * x[6] + x[7] * minus_x7) * _3;
		auto poly_3_unique = (poly_2_pos_shared + ((x[6] * x[7] + x[4] * x[5]) * _2) + x[3] * neg_two_thirds) * _3;
		// 4 5
		auto poly_4_pos_shared = x[10] * x[10] + x[11] * minus_x11;
		auto poly_4_unique = (poly_4_pos_shared + ((x[10] * minus_x11) * _2) + x[2] * x[2] + x[3] * minus_x3 + x[4] * neg_two_thirds) * _3;
		auto poly_5_unique = (poly_4_pos_shared + ((x[2] * x[3] + x[10] * x[11]) * _2) + x[5] * neg_two_thirds) * _3;
		// 6 7
		auto poly_6_pos_shared = (offset + x[2] * x[10] + x[3] * minus_x11);
		auto poly_6_neg_shared = (x[2] * x[11] + x[3] * x[10]);
		auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[6] * one_third) * _6;
		auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[7] * one_third) * _6;
		// 8 9
		auto poly_8_unique = (((x[0] * x[8] + x[1] * minus_x9)) + x[8] * one_third) * _6;
		auto poly_9_unique = (((x[0] * x[9] + x[1] * x[8])) + x[9] * one_third) * _6;
		// 10 11
		auto poly_10_unique = (((x[4] * x[6] + x[5] * minus_x7)) + x[10] * one_third) * _6;
		auto poly_11_unique = (((x[4] * x[7] + x[5] * x[6])) + x[11] * one_third) * _6;
		return ring.batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
	}

	/*
	// Incomplete functions - commented out until properly implemented
	cyclotomic_square(std::array<SmallElement, 12> & a) {
		// If both expand and reduce mult by 3 after mont before mask then <0, 1> output.
		// Inline everything??
		auto expanded = expand(a); // mem -> inline
		auto m1_res = fp12_flat_cyclotomic_sqr<1>(expanded);
		auto m1_res_reduced = batch_ele_reduce<1>(m1_res);
		auto m2_res = fp12_flat_cyclotomic_sqr<2>(a);
		auto m2_res_rns_reduced = reduce<false>(m2_res, m1_res_reduced);
		auto result = batch_ele_reduce<2, 3>(m2_res_rns_reduced);
	}

	cyclotomic_square(std::array<ExpandedElement, 12> & a) {
		// If both expand and reduce mult by 3 after mont before mask then <0, 1> output.
		fp12_flat_cyclotomic_sqr<1>(a.m1); // mem -> mem
		fp12_flat_cyclotomic_sqr<2>(a.m2); // mem -> inline
		reduce(a); // inline, mem -> mem
	}
	*/

	// While this is fine, it pushes the burden onto the cyclo square function
	// And although it allows two function calls, they will be unique code.
	/*
	cyclotomic_square(std::array<SmallElement, 12> & a) {
		// If both expand and reduce mult by 3 after mont before mask then <0, 1> output.
		fp12_flat_cyclotomic_sqr<1, true>(a); // mem -> expand inline -> mem
		fp12_flat_cyclotomic_sqr<2>(a); // mem -> reduce inline -> mem
	}


	cyclotomic_square(std::arrau<ExpandedElement, 12> & a) {
		fp12_flat_cyclotomic_sqr<1, false>(a);
	}
	*/
	
	inline StandardElement fp12_flat_frobenius_map1(const StandardElement &x, const Ring &ring) const {
		auto frob1_2_2 = x[2] * frob1_2;
		auto frob1_2_3 = x[3] * frob1_2;
		
		auto frob1_3_4 = x[4] * frob1_3;
		
		auto frob1_4_5 = x[5] * frob1_4;
		
		auto frob1_5_6 = x[6] * frob1_5;
		
		auto frob1_6_6 = x[6] * frob1_6;
		auto frob1_6_7 = x[7] * frob1_6;
		
		auto frob1_7_8 = x[8] * frob1_7;
		auto frob1_7_9 = x[9] * frob1_7;
		
		auto frob1_8_9 = x[9] * frob1_8;
		
		auto frob1_9_10 = x[10] * frob1_9;
		
		auto frob1_10_10 = x[10] * frob1_10;
		auto frob1_10_11 = x[11] * frob1_10;
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
		auto frob2_1_2 = x[2] * frob2_1;
		auto frob2_1_3 = x[3] * frob2_1;
		
		auto frob2_2_4 = x[4] * frob2_2;
		auto frob2_2_5 = x[5] * frob2_2;
		
		auto frob2_3_6 = x[6] * frob2_3;
		auto frob2_3_7 = x[7] * frob2_3;
		
		auto frob2_5_10 = x[10] * frob2_5;
		auto frob2_5_11 = x[11] * frob2_5;
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
		auto frob3_2_6 = x[6] * frob3_2;
		auto frob3_2_8 = x[8] * frob3_2;
		auto frob3_2_9 = x[9] * frob3_2;
		auto frob3_2_10 = x[10] * frob3_2;
		auto frob3_2_11 = x[11] * frob3_2;
		
		auto frob3_3_9 = x[9] * frob3_3;
		auto frob3_3_10 = x[10] * frob3_3;
		auto frob3_3_6 = x[6] * frob3_3;
		auto frob3_3_7 = x[7] * frob3_3;
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

	// inline StandardElement fp12_flat_conjugate(const StandardElement &x, const Ring &ring) const {
	// 	return StandardElement{x[0], x[1], x[2], x[3], x[4], x[5], ring.standard_negate(x[6]), ring.standard_negate(x[7]), ring.standard_negate(x[8]), ring.standard_negate(x[9]), ring.standard_negate(x[10]), ring.standard_negate(x[11])};
	// }

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
	auto poly_0_pos_shared = offset + ((x[2] * x[4] + x[3] * minus_x5 + x[6] * minus_x10 + x[7] * x[11]) * _2) + x[8] * minus_x8 + x[9] * x[9];
	auto poly_0_neg_shared = (x[2] * x[5] + x[3] * x[4] + x[6] * minus_x11 + x[7] * minus_x10 + x[8] * minus_x9) * _2;
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * _2;
	// 2 3
	auto poly_2_pos_shared = offset + ((x[8] * minus_x10 + x[9] * x[11]) * _2) + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_neg_shared = (x[4] * x[5] + x[8] * minus_x11 + x[9] * minus_x10) * _2;
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + ((x[0] * x[2] + x[1] * minus_x3) * _2) + x[6] * minus_x6 + x[7] * x[7];
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + (x[0] * x[3] + x[1] * x[2] + x[6] * minus_x7) * _2;
	// 4 5
	auto poly_4_pos_shared = offset + x[10] * minus_x10 + x[11] * x[11];
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
		auto poly_0_pos_shared = offset + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3;
		auto poly_0_neg_shared = x[2] * y[5] + x[3] * y[4] + x[4] * y[3] + x[5] * y[2];
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
		// 2 3
		auto poly_2_pos_shared = offset + x[4] * y[4] + x[5] * minus_y5;
		auto poly_2_neg_shared = x[4] * y[5] + x[5] * y[4];
		auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y[2] + x[1] * minus_y3 + x[2] * y[0] + x[3] * minus_y1;
		auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0];
		// 4 5
		auto poly_4_unique = x[0] * y[4] + x[1] * minus_y5 + x[2] * y[2] + x[3] * minus_y3 + x[4] * y[0] + x[5] * minus_y1;
		auto poly_5_unique = x[0] * y[5] + x[1] * y[4] + x[2] * y[3] + x[3] * y[2] + x[4] * y[1] + x[5] * y[0];
		// 6 7
		auto poly_6_pos_shared = offset + x[8] * minus_y4 + x[9] * y[5] + x[10] * minus_y2 + x[11] * y[3];
		auto poly_6_neg_shared = x[8] * minus_y5 + x[9] * minus_y4 + x[10] * minus_y3 + x[11] * minus_y2;
		auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[6] * minus_y0 + x[7] * y[1];
		auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[6] * minus_y1 + x[7] * minus_y0;
		// 8 9
		auto poly_8_pos_shared = offset + x[10] * minus_y4 + x[11] * y[5];
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
		auto poly_0_pos_shared = offset + x[2] * minus_x4 + x[3] * x[5];
		auto poly_0_neg_shared = x[2] * minus_x5 + x[3] * minus_x4;
		auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
		auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * _2;
		// 2 3
		auto poly_2_pos_shared = offset + x[4] * x[4] + x[5] * minus_x5;
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
		auto poly_0_pos_shared = offset + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3;
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
