/*
 * Copyright Supranational LLC
 * Copyright 2026 Peter Rindal
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include "blst/backend/backend.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fr.h"

using namespace blst;

#if defined(__x86_64__) || defined(__x86_64) || \
	(defined(_M_X64) && !defined(_M_ARM64EC))
# define CATBLST_AUDIT_X86_64 1
#else
# define CATBLST_AUDIT_X86_64 0
#endif

#if defined(BLST_ASM) && CATBLST_AUDIT_X86_64
extern "C" int blst_platform_cap;
#endif

namespace
{
#if defined(BLST_ASM) && CATBLST_AUDIT_X86_64
class PlatformCapabilityGuard
{
	int saved_;

public:
	explicit PlatformCapabilityGuard(int disabled_bits)
		: saved_(blst_platform_cap)
	{
		blst_platform_cap &= ~disabled_bits;
	}

	~PlatformCapabilityGuard()
	{
		blst_platform_cap = saved_;
	}
};
#endif

template<size_t N>
using Bytes = std::array<uint8_t, N>;

template<size_t N>
Bytes<N * sizeof(limb_t)> limbs_to_bytes(const limb_t (&in)[N])
{
	Bytes<N * sizeof(limb_t)> out{};
	for (size_t i = 0; i < N; ++i)
		for (size_t j = 0; j < sizeof(limb_t); ++j)
			out[i * sizeof(limb_t) + j] = uint8_t(in[i] >> (8 * j));
	return out;
}

template<size_t N>
void bytes_to_limbs(limb_t (&out)[N], const Bytes<N * sizeof(limb_t)>& in)
{
	for (size_t i = 0; i < N; ++i) {
		out[i] = 0;
		for (size_t j = 0; j < sizeof(limb_t); ++j)
			out[i] |= limb_t(in[i * sizeof(limb_t) + j]) << (8 * j);
	}
}

template<size_t N>
int compare(const Bytes<N>& a, const Bytes<N>& b)
{
	for (size_t i = N; i-- != 0;) {
		if (a[i] != b[i])
			return a[i] < b[i] ? -1 : 1;
	}
	return 0;
}

template<size_t N>
uint8_t add_raw(Bytes<N>& out, const Bytes<N>& a, const Bytes<N>& b)
{
	unsigned carry = 0;
	for (size_t i = 0; i < N; ++i) {
		const unsigned sum = unsigned(a[i]) + unsigned(b[i]) + carry;
		out[i] = uint8_t(sum);
		carry = sum >> 8;
	}
	return uint8_t(carry);
}

template<size_t N>
uint8_t sub_raw(Bytes<N>& out, const Bytes<N>& a, const Bytes<N>& b)
{
	unsigned borrow = 0;
	for (size_t i = 0; i < N; ++i) {
		const unsigned subtrahend = unsigned(b[i]) + borrow;
		out[i] = uint8_t(unsigned(a[i]) - subtrahend);
		borrow = unsigned(a[i]) < subtrahend;
	}
	return uint8_t(borrow);
}

template<size_t N>
Bytes<N> add_mod(const Bytes<N>& a, const Bytes<N>& b, const Bytes<N>& modulus)
{
	Bytes<N> sum{}, reduced{};
	const uint8_t carry = add_raw(sum, a, b);
	if (carry || compare(sum, modulus) >= 0) {
		(void)sub_raw(reduced, sum, modulus);
		return reduced;
	}
	return sum;
}

template<size_t N>
Bytes<N> sub_mod(const Bytes<N>& a, const Bytes<N>& b, const Bytes<N>& modulus)
{
	Bytes<N> difference{};
	if (sub_raw(difference, a, b) == 0)
		return difference;
	(void)add_raw(difference, difference, modulus);
	return difference;
}

template<size_t N>
Bytes<N> decrement(Bytes<N> value)
{
	for (size_t i = 0; i < N; ++i) {
		if (value[i]-- != 0)
			break;
	}
	return value;
}

template<size_t N>
Bytes<N> canonical_pattern(const Bytes<N>& modulus, uint8_t seed)
{
	Bytes<N> value{};
	for (size_t i = 0; i < N; ++i)
		value[i] = uint8_t(seed + 29 * i + (i >> 1));

	uint8_t top_mask = 0;
	for (uint8_t top = modulus[N - 1]; top != 0; top >>= 1)
		top_mask = uint8_t((top_mask << 1) | 1);
	value[N - 1] &= top_mask;
	while (compare(value, modulus) >= 0) {
		Bytes<N> reduced{};
		(void)sub_raw(reduced, value, modulus);
		value = reduced;
	}
	return value;
}

template<size_t N>
Bytes<N> lshift_mod(Bytes<N> value, size_t count, const Bytes<N>& modulus)
{
	while (count-- != 0)
		value = add_mod(value, value, modulus);
	return value;
}

template<size_t N>
Bytes<N> rshift_mod(Bytes<N> value, size_t count, const Bytes<N>& modulus)
{
	while (count-- != 0) {
		Bytes<N + 1> extended{};
		for (size_t i = 0; i < N; ++i)
			extended[i] = value[i];
		if ((value[0] & 1) != 0) {
			unsigned carry = 0;
			for (size_t i = 0; i < N; ++i) {
				const unsigned sum = unsigned(extended[i]) +
					unsigned(modulus[i]) + carry;
				extended[i] = uint8_t(sum);
				carry = sum >> 8;
			}
			extended[N] = uint8_t(carry);
		}
		for (size_t i = 0; i < N; ++i)
			value[i] = uint8_t((extended[i] >> 1) |
				(unsigned(extended[i + 1] & 1) << 7));
	}
	return value;
}

template<size_t DividendBytes, size_t DivisorBytes>
struct DivisionResult {
	Bytes<DividendBytes> quotient{};
	Bytes<DivisorBytes> remainder{};
};

template<size_t DividendBytes, size_t DivisorBytes>
DivisionResult<DividendBytes, DivisorBytes> divide_bits(
	const Bytes<DividendBytes>& dividend, const Bytes<DivisorBytes>& divisor)
{
	DivisionResult<DividendBytes, DivisorBytes> result;
	Bytes<DivisorBytes + 1> remainder{};
	Bytes<DivisorBytes + 1> extended_divisor{};
	for (size_t i = 0; i < DivisorBytes; ++i)
		extended_divisor[i] = divisor[i];

	for (size_t bit = DividendBytes * 8; bit-- != 0;) {
		unsigned carry = (dividend[bit / 8] >> (bit % 8)) & 1;
		for (size_t i = 0; i < remainder.size(); ++i) {
			const unsigned next = (unsigned(remainder[i]) << 1) | carry;
			remainder[i] = uint8_t(next);
			carry = next >> 8;
		}
		if (compare(remainder, extended_divisor) >= 0) {
			Bytes<DivisorBytes + 1> difference{};
			(void)sub_raw(difference, remainder, extended_divisor);
			remainder = difference;
			result.quotient[bit / 8] |= uint8_t(1u << (bit % 8));
		}
	}
	for (size_t i = 0; i < DivisorBytes; ++i)
		result.remainder[i] = remainder[i];
	return result;
}

template<size_t N>
uint64_t low_u64(const Bytes<N>& value)
{
	static_assert(N >= sizeof(uint64_t));
	uint64_t out = 0;
	for (size_t i = 0; i < sizeof(out); ++i)
		out |= uint64_t(value[i]) << (8 * i);
	return out;
}

template<size_t A, size_t B>
Bytes<A + B> multiply_wide(const Bytes<A>& a, const Bytes<B>& b)
{
	Bytes<A + B> product{};
	for (size_t i = 0; i < A; ++i) {
		unsigned carry = 0;
		for (size_t j = 0; j < B; ++j) {
			const unsigned term = unsigned(product[i + j]) +
				unsigned(a[i]) * unsigned(b[j]) + carry;
			product[i + j] = uint8_t(term);
			carry = term >> 8;
		}
		for (size_t k = i + B; carry != 0; ++k) {
			if (k >= product.size()) {
				ADD_FAILURE() << "byte multiplication overflowed its exact width";
				return product;
			}
			const unsigned term = unsigned(product[k]) + carry;
			product[k] = uint8_t(term);
			carry = term >> 8;
		}
	}
	return product;
}

template<size_t DividendBytes, size_t ModulusBytes>
Bytes<ModulusBytes> reduce_mod(const Bytes<DividendBytes>& value,
	const Bytes<ModulusBytes>& modulus)
{
	return divide_bits(value, modulus).remainder;
}

template<size_t N>
Bytes<N> multiply_mod(const Bytes<N>& a, const Bytes<N>& b,
	const Bytes<N>& modulus)
{
	return reduce_mod(multiply_wide(a, b), modulus);
}

template<size_t N>
Bytes<N> radix_mod(const Bytes<N>& modulus)
{
	Bytes<N + 1> radix{};
	radix[N] = 1;
	return reduce_mod(radix, modulus);
}

template<size_t N>
Bytes<N> left_align(const Bytes<N>& value)
{
	unsigned top_bits = 0;
	for (uint8_t top = value[N - 1]; top != 0; top >>= 1)
		++top_bits;
	EXPECT_NE(top_bits, 0u);
	const unsigned shift = 8 - top_bits;
	Bytes<N> out{};
	unsigned carry = 0;
	for (size_t i = 0; i < N; ++i) {
		const unsigned next = (unsigned(value[i]) << shift) | carry;
		out[i] = uint8_t(next);
		carry = next >> 8;
	}
	EXPECT_EQ(carry, 0u);
	return out;
}

template<size_t N>
Bytes<N> rshift_one(Bytes<N> value)
{
	unsigned carry = 0;
	for (size_t i = N; i-- != 0;) {
		const unsigned next = value[i];
		value[i] = uint8_t((next >> 1) | (carry << 7));
		carry = next & 1;
	}
	return value;
}

template<size_t N>
Bytes<N> pow_mod(Bytes<N> base, const Bytes<N>& exponent,
	const Bytes<N>& modulus)
{
	Bytes<N> result{};
	result[0] = 1;
	for (size_t byte_index = 0; byte_index < N; ++byte_index) {
		for (unsigned bit = 0; bit < 8; ++bit) {
			if ((exponent[byte_index] & (1u << bit)) != 0)
				result = multiply_mod(result, base, modulus);
			base = multiply_mod(base, base, modulus);
		}
	}
	return result;
}

template<size_t Limbs, class Inverse>
void expect_almost_inverse_provider(const char* provider,
	const limb_t (&modulus_limbs)[Limbs], Inverse&& inverse)
{
	SCOPED_TRACE(provider);
	constexpr size_t N = Limbs * sizeof(limb_t);
	const Bytes<N> modulus = limbs_to_bytes(modulus_limbs);
	const Bytes<N> modulus_x = left_align(modulus);
	const Bytes<N> radix = radix_mod(modulus);
	const Bytes<N> radix_squared = multiply_mod(radix, radix, modulus);
	Bytes<N> one{};
	one[0] = 1;
	const std::array<Bytes<N>, 4> inputs = {
		Bytes<N>{}, one, decrement(modulus), canonical_pattern(modulus, 0x3d)
	};

	Bytes<2 * N> bound{};
	for (size_t i = 0; i < N; ++i)
		bound[N + i] = modulus_x[i];

	limb_t modx_limbs[Limbs];
	bytes_to_limbs(modx_limbs, modulus_x);
	for (const Bytes<N>& input : inputs) {
		limb_t input_limbs[Limbs], output_limbs[2 * Limbs];
		bytes_to_limbs(input_limbs, input);
		inverse(output_limbs, input_limbs, modulus_limbs, modx_limbs);
		const Bytes<2 * N> output = limbs_to_bytes(output_limbs);
		EXPECT_LT(compare(output, bound), 0);
		if (compare(input, Bytes<N>{}) == 0) {
			EXPECT_EQ(output, Bytes<2 * N>{});
			continue;
		}
		const Bytes<N> reduced = reduce_mod(output, modulus);
		EXPECT_EQ(multiply_mod(reduced, input, modulus), radix_squared);
	}
}

bool expected_is_square(const Bytes<48>& input, const Bytes<48>& modulus)
{
	const Bytes<48> zero{};
	if (input == zero)
		return true;
	Bytes<48> one{};
	one[0] = 1;
	if (input == one)
		return true;
	const Bytes<48> exponent = rshift_one(decrement(modulus));
	const Bytes<48> symbol = pow_mod(input, exponent, modulus);
	const Bytes<48> minus_one = decrement(modulus);
	EXPECT_TRUE(symbol == one || symbol == minus_one);
	return symbol == one;
}

template<class IsSquare>
void expect_is_square_provider(const char* provider, IsSquare&& is_square)
{
	SCOPED_TRACE(provider);
	const Bytes<48> modulus = limbs_to_bytes(get_BLST_P().l);
	Bytes<48> one{};
	one[0] = 1;
	const Bytes<48> pattern = canonical_pattern(modulus, 0x71);
	const Bytes<48> square = multiply_mod(pattern, pattern, modulus);
	const std::array<Bytes<48>, 5> inputs = {
		Bytes<48>{}, one, decrement(modulus), pattern, square
	};
	for (const Bytes<48>& input : inputs) {
		Fp input_limbs{};
		bytes_to_limbs(input_limbs.l, input);
		EXPECT_EQ(is_square(input_limbs), expected_is_square(input, modulus));
	}
}

template<size_t N>
void expect_montgomery_product(const Bytes<N>& out, const Bytes<N>& a,
	const Bytes<N>& b, const Bytes<N>& modulus)
{
	EXPECT_LT(compare(out, modulus), 0);
	EXPECT_EQ(multiply_mod(out, radix_mod(modulus), modulus),
		multiply_mod(a, b, modulus));
}

template<size_t N, size_t WideN>
void expect_montgomery_reduction(const Bytes<N>& out,
	const Bytes<WideN>& input, const Bytes<N>& modulus)
{
	EXPECT_LT(compare(out, modulus), 0);
	EXPECT_EQ(multiply_mod(out, radix_mod(modulus), modulus),
		reduce_mod(input, modulus));
}

template<size_t N>
void expect_montgomery_value(const Bytes<N>& out, const Bytes<N>& value,
	const Bytes<N>& modulus)
{
	EXPECT_LT(compare(out, modulus), 0);
	EXPECT_EQ(multiply_mod(out, radix_mod(modulus), modulus), value);
}

template<size_t N>
std::array<Bytes<N>, 4> montgomery_values(const Bytes<N>& modulus)
{
	Bytes<N> one{};
	one[0] = 1;
	return {
		Bytes<N>{}, one, decrement(modulus),
		canonical_pattern(modulus, 0x6d)
	};
}

template<size_t Limbs, class Mul, class Sqr, class Redc, class From>
void expect_montgomery_provider(const char* provider,
	const limb_t (&modulus_limbs)[Limbs], limb_t n0,
	bool sparse_second_operand, Mul&& mul, Sqr&& sqr, Redc&& redc, From&& from)
{
	SCOPED_TRACE(provider);
	constexpr size_t N = Limbs * sizeof(limb_t);
	const Bytes<N> modulus = limbs_to_bytes(modulus_limbs);
	const auto values = montgomery_values(modulus);
	Bytes<N> full_width{};
	full_width.fill(0xff);

	for (size_t i = 0; i < values.size(); ++i) {
		const Bytes<N>& a_bytes = values[i];
		const Bytes<N>& b_bytes = sparse_second_operand && i + 1 == values.size()
			? full_width : values[(i + 1) % values.size()];
		limb_t a[Limbs], b[Limbs], out[Limbs], alias[Limbs];
		bytes_to_limbs(a, a_bytes);
		bytes_to_limbs(b, b_bytes);

		mul(out, a, b, modulus_limbs, n0);
		expect_montgomery_product(limbs_to_bytes(out), a_bytes, b_bytes,
			modulus);
		std::memcpy(alias, a, sizeof(alias));
		mul(alias, alias, b, modulus_limbs, n0);
		expect_montgomery_product(limbs_to_bytes(alias), a_bytes, b_bytes,
			modulus);
		std::memcpy(alias, b, sizeof(alias));
		mul(alias, a, alias, modulus_limbs, n0);
		expect_montgomery_product(limbs_to_bytes(alias), a_bytes, b_bytes,
			modulus);

		sqr(out, a, modulus_limbs, n0);
		expect_montgomery_product(limbs_to_bytes(out), a_bytes, a_bytes,
			modulus);
		std::memcpy(alias, a, sizeof(alias));
		sqr(alias, alias, modulus_limbs, n0);
		expect_montgomery_product(limbs_to_bytes(alias), a_bytes, a_bytes,
			modulus);

		from(out, a, modulus_limbs, n0);
		expect_montgomery_reduction(limbs_to_bytes(out), a_bytes, modulus);
		std::memcpy(alias, a, sizeof(alias));
		from(alias, alias, modulus_limbs, n0);
		expect_montgomery_reduction(limbs_to_bytes(alias), a_bytes, modulus);
	}
	{
		limb_t a[Limbs], out[Limbs], alias[Limbs];
		bytes_to_limbs(a, full_width);
		from(out, a, modulus_limbs, n0);
		expect_montgomery_reduction(limbs_to_bytes(out), full_width, modulus);
		std::memcpy(alias, a, sizeof(alias));
		from(alias, alias, modulus_limbs, n0);
		expect_montgomery_reduction(limbs_to_bytes(alias), full_width, modulus);
	}

	for (size_t i = 0; i < values.size(); ++i) {
		Bytes<2 * N> input{};
		const Bytes<N>& upper = values[i];
		for (size_t j = 0; j < N; ++j) {
			input[j] = uint8_t(0xa7 + 31 * j + i);
			input[N + j] = upper[j];
		}
		limb_t wide[2 * Limbs], out[Limbs];
		bytes_to_limbs(wide, input);
		redc(out, wide, modulus_limbs, n0);
		expect_montgomery_reduction(limbs_to_bytes(out), input, modulus);

		limb_t alias[2 * Limbs];
		std::memcpy(alias, wide, sizeof(alias));
		redc(alias, alias, modulus_limbs, n0);
		limb_t alias_out[Limbs];
		std::memcpy(alias_out, alias, sizeof(alias_out));
		expect_montgomery_reduction(limbs_to_bytes(alias_out), input, modulus);
	}
}

template<class Mul, class Sqr>
void expect_fp2_montgomery_provider(const char* provider, Mul&& mul, Sqr&& sqr)
{
	SCOPED_TRACE(provider);
	const Bytes<48> modulus = limbs_to_bytes(get_BLST_P().l);
	const auto values = montgomery_values(modulus);
	for (size_t i = 0; i < 3; ++i) {
		const Bytes<48>& a0 = values[i + 1];
		const Bytes<48>& a1 = values[(i + 2) % values.size()];
		const Bytes<48>& b0 = values[(i + 3) % values.size()];
		const Bytes<48>& b1 = values[i];
		Fp2 a{}, b{}, out{}, alias{};
		bytes_to_limbs(a[0].l, a0);
		bytes_to_limbs(a[1].l, a1);
		bytes_to_limbs(b[0].l, b0);
		bytes_to_limbs(b[1].l, b1);

		const Bytes<48> real = sub_mod(multiply_mod(a0, b0, modulus),
			multiply_mod(a1, b1, modulus), modulus);
		const Bytes<48> imaginary = add_mod(multiply_mod(a0, b1, modulus),
			multiply_mod(a1, b0, modulus), modulus);
		mul(out, a, b);
		expect_montgomery_value(limbs_to_bytes(out[0].l), real, modulus);
		expect_montgomery_value(limbs_to_bytes(out[1].l), imaginary, modulus);
		alias = a;
		mul(alias, alias, b);
		expect_montgomery_value(limbs_to_bytes(alias[0].l), real, modulus);
		expect_montgomery_value(limbs_to_bytes(alias[1].l), imaginary, modulus);
		alias = b;
		mul(alias, a, alias);
		expect_montgomery_value(limbs_to_bytes(alias[0].l), real, modulus);
		expect_montgomery_value(limbs_to_bytes(alias[1].l), imaginary, modulus);

		const Bytes<48> square_real = sub_mod(multiply_mod(a0, a0, modulus),
			multiply_mod(a1, a1, modulus), modulus);
		const Bytes<48> square_imaginary = add_mod(
			multiply_mod(a0, a1, modulus), multiply_mod(a0, a1, modulus),
			modulus);
		sqr(out, a);
		expect_montgomery_value(limbs_to_bytes(out[0].l), square_real, modulus);
		expect_montgomery_value(limbs_to_bytes(out[1].l), square_imaginary,
			modulus);
		alias = a;
		sqr(alias, alias);
		expect_montgomery_value(limbs_to_bytes(alias[0].l), square_real,
			modulus);
		expect_montgomery_value(limbs_to_bytes(alias[1].l), square_imaginary,
			modulus);
	}
}

template<class Mul, class Sqr>
void expect_fp2_wide_provider(const char* provider, Mul&& mul, Sqr&& sqr)
{
	SCOPED_TRACE(provider);
	const Bytes<48> modulus = limbs_to_bytes(get_BLST_P().l);
	const auto values = montgomery_values(modulus);
	for (size_t i = 0; i < 3; ++i) {
		const Bytes<48>& a0 = values[i + 1];
		const Bytes<48>& a1 = values[(i + 2) % values.size()];
		const Bytes<48>& b0 = values[(i + 3) % values.size()];
		const Bytes<48>& b1 = values[i];
		Fp2 a{}, b{};
		UFp2 out{};
		bytes_to_limbs(a[0].l, a0);
		bytes_to_limbs(a[1].l, a1);
		bytes_to_limbs(b[0].l, b0);
		bytes_to_limbs(b[1].l, b1);

		const Bytes<48> real = sub_mod(multiply_mod(a0, b0, modulus),
			multiply_mod(a1, b1, modulus), modulus);
		const Bytes<48> imaginary = add_mod(multiply_mod(a0, b1, modulus),
			multiply_mod(a1, b0, modulus), modulus);
		mul(out, a, b);
		EXPECT_EQ(reduce_mod(limbs_to_bytes(out[0].data), modulus), real);
		EXPECT_EQ(reduce_mod(limbs_to_bytes(out[1].data), modulus), imaginary);

		const Bytes<48> square_real = sub_mod(multiply_mod(a0, a0, modulus),
			multiply_mod(a1, a1, modulus), modulus);
		const Bytes<48> square_imaginary = add_mod(
			multiply_mod(a0, a1, modulus), multiply_mod(a0, a1, modulus),
			modulus);
		sqr(out, a);
		EXPECT_EQ(reduce_mod(limbs_to_bytes(out[0].data), modulus), square_real);
		EXPECT_EQ(reduce_mod(limbs_to_bytes(out[1].data), modulus),
			square_imaginary);
	}
}

template<class SqrNMul>
void expect_repeated_square_provider(const char* provider, SqrNMul&& operation)
{
	SCOPED_TRACE(provider);
	const Bytes<48> modulus = limbs_to_bytes(get_BLST_P().l);
	const auto values = montgomery_values(modulus);
	const Bytes<48> a_bytes = values[3];
	const Bytes<48> b_bytes = values[2];
	for (size_t count : { size_t(1), size_t(2), size_t(5) }) {
		Fp a{}, b{}, out{}, alias{};
		bytes_to_limbs(a.l, a_bytes);
		bytes_to_limbs(b.l, b_bytes);
		operation(out, a, count, b);

		Bytes<48> power = a_bytes;
		Bytes<48> radix_power = radix_mod(modulus);
		for (size_t i = 0; i < count; ++i) {
			power = multiply_mod(power, power, modulus);
			radix_power = multiply_mod(radix_power, radix_power, modulus);
		}
		const Bytes<48> expected = multiply_mod(power, b_bytes, modulus);
		const Bytes<48> out_bytes = limbs_to_bytes(out.l);
		EXPECT_LT(compare(out_bytes, modulus), 0);
		EXPECT_EQ(multiply_mod(out_bytes, radix_power, modulus), expected);

		alias = a;
		operation(alias, alias, count, b);
		const Bytes<48> alias_bytes = limbs_to_bytes(alias.l);
		EXPECT_LT(compare(alias_bytes, modulus), 0);
		EXPECT_EQ(multiply_mod(alias_bytes, radix_power, modulus), expected);
	}
}

template<size_t Limbs>
void expect_single_width_modular_arithmetic(const limb_t (&modulus_limbs)[Limbs])
{
	constexpr size_t N = Limbs * sizeof(limb_t);
	const Bytes<N> modulus = limbs_to_bytes(modulus_limbs);
	const Bytes<N> zero{};
	Bytes<N> one{};
	one[0] = 1;
	const Bytes<N> modulus_minus_one = decrement(modulus);
	const std::array<Bytes<N>, 5> values = {
		zero, one, modulus_minus_one, canonical_pattern(modulus, 0x37),
		canonical_pattern(modulus, 0xa1)
	};

	for (const Bytes<N>& a_bytes : values) {
		for (const Bytes<N>& b_bytes : values) {
			limb_t a[Limbs], b[Limbs], out[Limbs], alias[Limbs];
			bytes_to_limbs(a, a_bytes);
			bytes_to_limbs(b, b_bytes);
			const Bytes<N> expected_add = add_mod(a_bytes, b_bytes, modulus);
			const Bytes<N> expected_sub = sub_mod(a_bytes, b_bytes, modulus);

			if constexpr (Limbs == NLIMBS(256)) {
				add_mod_256(out, a, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(out), expected_add);
				std::memcpy(alias, a, sizeof(alias));
				add_mod_256(alias, alias, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_add);
				std::memcpy(alias, b, sizeof(alias));
				add_mod_256(alias, a, alias, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_add);
				sub_mod_256(out, a, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(out), expected_sub);
				std::memcpy(alias, a, sizeof(alias));
				sub_mod_256(alias, alias, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_sub);
				std::memcpy(alias, b, sizeof(alias));
				sub_mod_256(alias, a, alias, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_sub);
			} else {
				add_mod_384(out, a, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(out), expected_add);
				std::memcpy(alias, a, sizeof(alias));
				add_mod_384(alias, alias, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_add);
				std::memcpy(alias, b, sizeof(alias));
				add_mod_384(alias, a, alias, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_add);
				sub_mod_384(out, a, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(out), expected_sub);
				std::memcpy(alias, a, sizeof(alias));
				sub_mod_384(alias, alias, b, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_sub);
				std::memcpy(alias, b, sizeof(alias));
				sub_mod_384(alias, a, alias, modulus_limbs);
				EXPECT_EQ(limbs_to_bytes(alias), expected_sub);
			}
		}
	}

	for (size_t count : { size_t(1), size_t(3), size_t(64) }) {
		limb_t a[Limbs], out[Limbs], alias[Limbs];
		bytes_to_limbs(a, modulus_minus_one);
		if constexpr (Limbs == NLIMBS(256)) {
			lshift_mod_256(out, a, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(out), lshift_mod(modulus_minus_one, count, modulus));
			std::memcpy(alias, a, sizeof(alias));
			lshift_mod_256(alias, alias, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(alias), lshift_mod(modulus_minus_one, count, modulus));
			rshift_mod_256(out, a, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(out), rshift_mod(modulus_minus_one, count, modulus));
			std::memcpy(alias, a, sizeof(alias));
			rshift_mod_256(alias, alias, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(alias), rshift_mod(modulus_minus_one, count, modulus));
		} else {
			lshift_mod_384(out, a, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(out), lshift_mod(modulus_minus_one, count, modulus));
			std::memcpy(alias, a, sizeof(alias));
			lshift_mod_384(alias, alias, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(alias), lshift_mod(modulus_minus_one, count, modulus));
			rshift_mod_384(out, a, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(out), rshift_mod(modulus_minus_one, count, modulus));
			std::memcpy(alias, a, sizeof(alias));
			rshift_mod_384(alias, alias, count, modulus_limbs);
			EXPECT_EQ(limbs_to_bytes(alias), rshift_mod(modulus_minus_one, count, modulus));
		}
	}
}
} // namespace

TEST(AssemblyArithmeticAudit, ZeroLengthSelectedPrefetchIsNoOp)
{
	vec_prefetch(nullptr, 0);
}

#ifdef BLST_ASM
TEST(AssemblyArithmeticAudit, SingleWidthOperationsMatchByteOracle)
{
	expect_single_width_modular_arithmetic(get_BLST_r().l);
	expect_single_width_modular_arithmetic(get_BLST_P().l);
}

TEST(AssemblyArithmeticAudit, SelectedMontgomeryProvidersMatchByteOracle)
{
	expect_montgomery_provider("selected sparse 256", get_BLST_r().l,
		BLST_r0, true,
		[](limb_t* ret, const limb_t* a, const limb_t* b,
			const limb_t* p, limb_t n0) {
			mul_mont_sparse_v256(ret, a, b, p, n0);
		},
		[](limb_t* ret, const limb_t* a, const limb_t* p, limb_t n0) {
			sqr_mont_sparse_v256(ret, a, p, n0);
		},
		[](limb_t* ret, const limb_t* a, const limb_t* p, limb_t n0) {
			redc_mont_v256(ret, a, p, n0);
		},
		[](limb_t* ret, const limb_t* a, const limb_t* p, limb_t n0) {
			from_mont_v256(ret, a, p, n0);
		});

	expect_montgomery_provider("selected 384", get_BLST_P().l,
		BLST_p0, false,
		[](limb_t* ret, const limb_t* a, const limb_t* b,
			const limb_t* p, limb_t n0) {
			mul_mont_v384(ret, a, b, p, n0);
		},
		[](limb_t* ret, const limb_t* a, const limb_t* p, limb_t n0) {
			sqr_mont_v384(ret, a, p, n0);
		},
		[](limb_t* ret, const limb_t* a, const limb_t* p, limb_t n0) {
			redc_mont_v384(ret, a, p, n0);
		},
		[](limb_t* ret, const limb_t* a, const limb_t* p, limb_t n0) {
			from_mont_v384(ret, a, p, n0);
		});
}

TEST(AssemblyArithmeticAudit, SelectedAlmostInverseMatchesByteOracle)
{
	expect_almost_inverse_provider("selected 256 almost-inverse",
		get_BLST_r().l,
		[](limb_t* ret, const limb_t* inp, const limb_t* mod,
			const limb_t* modx) {
			ct_inverse_mod_v256(ret, inp, mod, modx);
		});
	expect_almost_inverse_provider("selected 384 almost-inverse",
		get_BLST_P().l,
		[](limb_t* ret, const limb_t* inp, const limb_t* mod,
			const limb_t* modx) {
			ct_inverse_mod_v384(ret, inp, mod, modx);
		});
}

TEST(AssemblyArithmeticAudit, PortableAlmostInverseMatchesByteOracle)
{
	expect_almost_inverse_provider("portable 256 almost-inverse",
		get_BLST_r().l,
		[](limb_t* ret, const limb_t* inp, const limb_t* mod,
			const limb_t* modx) {
			cc_ct_inverse_mod_n<NLIMBS(256)>(ret, inp, mod, modx);
		});
	expect_almost_inverse_provider("portable 384 almost-inverse",
		get_BLST_P().l,
		[](limb_t* ret, const limb_t* inp, const limb_t* mod,
			const limb_t* modx) {
			cc_ct_inverse_mod_n<NLIMBS(384)>(ret, inp, mod, modx);
		});
}

#if CATBLST_AUDIT_X86_64 && \
	(defined(_MSC_VER) || defined(BLST_ADX) || defined(BLST_PORTABLE))
TEST(AssemblyArithmeticAudit, AdxAlmostInverseMatchesByteOracle)
{
	expect_almost_inverse_provider("ADX 384 almost-inverse", get_BLST_P().l,
		ctx_inverse_mod_384);
}
#endif

#if (!CATBLST_AUDIT_X86_64) || defined(_MSC_VER) || \
	!defined(BLST_ADX) || defined(BLST_PORTABLE)
TEST(AssemblyArithmeticAudit, GenericAlmostInverseMatchesByteOracle)
{
#if CATBLST_AUDIT_X86_64
	PlatformCapabilityGuard force_generic(1);
#endif
	expect_almost_inverse_provider("generic 384 almost-inverse", get_BLST_P().l,
		ct_inverse_mod_384);
}
#endif

TEST(AssemblyArithmeticAudit, IsSquareMatchesEulerCriterion)
{
	expect_is_square_provider("assembly quadratic-residue test",
		[](const Fp& input) {
			return ct_is_square_mod_384(input.l, get_BLST_P().l) != 0;
		});
	expect_is_square_provider("portable quadratic-residue test",
		[](const Fp& input) {
			return cc_ct_is_sqr_mod_n<NLIMBS(384)>(input.l,
				get_BLST_P().l) != 0;
		});
}

#if CATBLST_AUDIT_X86_64 && \
	(defined(_MSC_VER) || defined(BLST_ADX) || defined(BLST_PORTABLE))
TEST(AssemblyArithmeticAudit, AdxMontgomeryProvidersMatchByteOracle)
{
	expect_montgomery_provider("ADX sparse 256", get_BLST_r().l,
		BLST_r0, true, mulx_mont_sparse_256, sqrx_mont_sparse_256,
		redcx_mont_256, fromx_mont_256);
	expect_montgomery_provider("ADX 384", get_BLST_P().l,
		BLST_p0, false, mulx_mont_384, sqrx_mont_384,
		redcx_mont_384, fromx_mont_384);
}
#endif

#if (!CATBLST_AUDIT_X86_64) || defined(_MSC_VER) || \
	!defined(BLST_ADX) || defined(BLST_PORTABLE)
TEST(AssemblyArithmeticAudit, GenericMontgomeryProvidersMatchByteOracle)
{
#if CATBLST_AUDIT_X86_64
	PlatformCapabilityGuard force_generic(1);
#endif
	expect_montgomery_provider("generic sparse 256", get_BLST_r().l,
		BLST_r0, true, mul_mont_sparse_256, sqr_mont_sparse_256,
		redc_mont_256, from_mont_256);
	expect_montgomery_provider("generic 384", get_BLST_P().l,
		BLST_p0, false, mul_mont_384, sqr_mont_384,
		redc_mont_384, from_mont_384);
}
#endif

TEST(AssemblyArithmeticAudit, SelectedFp2MontgomeryProvidersMatchByteOracle)
{
	expect_fp2_montgomery_provider("selected Fp2",
		[](Fp2& ret, const Fp2& a, const Fp2& b) {
			Fp2::mul(ret, a, b);
		},
		[](Fp2& ret, const Fp2& a) {
			Fp2::sqr(ret, a);
		});
	expect_fp2_wide_provider("selected wide Fp2",
		[](UFp2& ret, const Fp2& a, const Fp2& b) {
			mul_v382x(ret, a, b, get_BLST_P());
		},
		[](UFp2& ret, const Fp2& a) {
			sqr_v382x(ret, a, get_BLST_P());
		});
}

#if CATBLST_AUDIT_X86_64 && \
	(defined(_MSC_VER) || defined(BLST_ADX) || defined(BLST_PORTABLE))
TEST(AssemblyArithmeticAudit, AdxFp2MontgomeryProvidersMatchByteOracle)
{
	expect_fp2_montgomery_provider("ADX Fp2",
		[](Fp2& ret, const Fp2& a, const Fp2& b) {
			mulx_mont_384x(&ret, &a, &b, &get_BLST_P(), BLST_p0);
		},
		[](Fp2& ret, const Fp2& a) {
			sqrx_mont_384x(&ret, &a, &get_BLST_P(), BLST_p0);
		});
	expect_fp2_wide_provider("ADX wide Fp2",
		[](UFp2& ret, const Fp2& a, const Fp2& b) {
			mulx_382x(&ret, &a, &b, &get_BLST_P());
		},
		[](UFp2& ret, const Fp2& a) {
			sqrx_382x(&ret, &a, &get_BLST_P());
		});
}
#endif

#if (!CATBLST_AUDIT_X86_64) || defined(_MSC_VER) || \
	!defined(BLST_ADX) || defined(BLST_PORTABLE)
TEST(AssemblyArithmeticAudit, GenericFp2MontgomeryProvidersMatchByteOracle)
{
#if CATBLST_AUDIT_X86_64
	PlatformCapabilityGuard force_generic(1);
#endif
	expect_fp2_montgomery_provider("generic Fp2",
		[](Fp2& ret, const Fp2& a, const Fp2& b) {
			mul_mont_384x(&ret, &a, &b, &get_BLST_P(), BLST_p0);
		},
		[](Fp2& ret, const Fp2& a) {
			sqr_mont_384x(&ret, &a, &get_BLST_P(), BLST_p0);
		});
	expect_fp2_wide_provider("generic wide Fp2",
		[](UFp2& ret, const Fp2& a, const Fp2& b) {
			mul_382x(&ret, &a, &b, &get_BLST_P());
		},
		[](UFp2& ret, const Fp2& a) {
			sqr_382x(&ret, &a, &get_BLST_P());
		});
}
#endif

TEST(AssemblyArithmeticAudit, SelectedRepeatedSquareMatchesByteOracle)
{
	expect_repeated_square_provider("selected repeated square",
		[](Fp& ret, const Fp& a, size_t count, const Fp& b) {
			Fp::sqr_n_mul(ret, a, count, b);
		});

	const Bytes<48> modulus = limbs_to_bytes(get_BLST_P().l);
	const auto values = montgomery_values(modulus);
	Fp a{}, b{}, out{};
	bytes_to_limbs(a.l, values[3]);
	bytes_to_limbs(b.l, values[2]);
	Fp::sqr_n_mul(out, a, 0, b);
	expect_montgomery_product(limbs_to_bytes(out.l), values[3], values[2],
		modulus);
	out = b;
	Fp::sqr_n_mul(out, a, 2, out);
	Bytes<48> power = multiply_mod(values[3], values[3], modulus);
	power = multiply_mod(power, power, modulus);
	Bytes<48> radix_power = multiply_mod(radix_mod(modulus),
		radix_mod(modulus), modulus);
	radix_power = multiply_mod(radix_power, radix_power, modulus);
	EXPECT_EQ(multiply_mod(limbs_to_bytes(out.l), radix_power, modulus),
		multiply_mod(power, values[2], modulus));
}

#if CATBLST_AUDIT_X86_64 && \
	(defined(_MSC_VER) || defined(BLST_ADX) || defined(BLST_PORTABLE))
TEST(AssemblyArithmeticAudit, AdxRepeatedSquareMatchesByteOracle)
{
	expect_repeated_square_provider("ADX repeated square",
		[](Fp& ret, const Fp& a, size_t count, const Fp& b) {
			sqrx_n_mul_mont_383(&ret, &a, count, &get_BLST_P(), BLST_p0, &b);
		});
}
#endif

#if (!CATBLST_AUDIT_X86_64) || defined(_MSC_VER) || \
	!defined(BLST_ADX) || defined(BLST_PORTABLE)
TEST(AssemblyArithmeticAudit, GenericRepeatedSquareMatchesByteOracle)
{
#if CATBLST_AUDIT_X86_64
	PlatformCapabilityGuard force_generic(1);
#endif
	expect_repeated_square_provider("generic repeated square",
		[](Fp& ret, const Fp& a, size_t count, const Fp& b) {
			sqr_n_mul_mont_383(&ret, &a, count, &get_BLST_P(), BLST_p0, &b);
		});
}
#endif

TEST(AssemblyArithmeticAudit, ScalarChecksCoverCanonicalBoundariesAndAliases)
{
	const Bytes<32> modulus = limbs_to_bytes(get_BLST_r().l);
	const Bytes<32> zero{};
	Bytes<32> one{};
	one[0] = 1;
	const Bytes<32> modulus_minus_one = decrement(modulus);
	Bytes<32> modulus_plus_one = modulus;
	modulus_plus_one[0]++;

	EXPECT_EQ(check_mod_256(zero.data(), get_BLST_r().l), 0u);
	EXPECT_EQ(check_mod_256(one.data(), get_BLST_r().l), 1u);
	EXPECT_EQ(check_mod_256(modulus_minus_one.data(), get_BLST_r().l), 1u);
	EXPECT_EQ(check_mod_256(modulus.data(), get_BLST_r().l), 0u);
	EXPECT_EQ(check_mod_256(modulus_plus_one.data(), get_BLST_r().l), 0u);

	Bytes<32> result = modulus_minus_one;
	EXPECT_EQ(add_n_check_mod_256(result.data(), result.data(), one.data(),
		get_BLST_r().l), 0u);
	EXPECT_EQ(result, zero);
	result = zero;
	EXPECT_EQ(sub_n_check_mod_256(result.data(), one.data(), result.data(),
		get_BLST_r().l), 1u);
	EXPECT_EQ(result, one);
}

TEST(AssemblyArithmeticAudit, DoubleWidthOperationsMatchByteOracle)
{
	Bytes<96> modulus{};
	const Bytes<48> field_modulus = limbs_to_bytes(get_BLST_P().l);
	for (size_t i = 0; i < field_modulus.size(); ++i)
		modulus[48 + i] = field_modulus[i];
	const Bytes<96> zero{};
	Bytes<96> one{};
	one[0] = 1;
	const Bytes<96> modulus_minus_one = decrement(modulus);
	Bytes<96> pattern{};
	const Bytes<48> upper_pattern = canonical_pattern(field_modulus, 0x5b);
	for (size_t i = 0; i < 48; ++i) {
		pattern[i] = uint8_t(0xc3 + 17 * i);
		pattern[48 + i] = upper_pattern[i];
	}
	const std::array<Bytes<96>, 4> values = {
		zero, one, modulus_minus_one, pattern
	};

	for (const Bytes<96>& a_bytes : values) {
		for (const Bytes<96>& b_bytes : values) {
			limb_t a[12], b[12], out[12], alias[12];
			bytes_to_limbs(a, a_bytes);
			bytes_to_limbs(b, b_bytes);
			const Bytes<96> expected_add = add_mod(a_bytes, b_bytes, modulus);
			const Bytes<96> expected_sub = sub_mod(a_bytes, b_bytes, modulus);

			add_mod_384x384(out, a, b, get_BLST_P().l);
			EXPECT_EQ(limbs_to_bytes(out), expected_add);
			std::memcpy(alias, a, sizeof(alias));
			add_mod_384x384(alias, alias, b, get_BLST_P().l);
			EXPECT_EQ(limbs_to_bytes(alias), expected_add);
			std::memcpy(alias, b, sizeof(alias));
			add_mod_384x384(alias, a, alias, get_BLST_P().l);
			EXPECT_EQ(limbs_to_bytes(alias), expected_add);
			sub_mod_384x384(out, a, b, get_BLST_P().l);
			EXPECT_EQ(limbs_to_bytes(out), expected_sub);
			std::memcpy(alias, a, sizeof(alias));
			sub_mod_384x384(alias, alias, b, get_BLST_P().l);
			EXPECT_EQ(limbs_to_bytes(alias), expected_sub);
			std::memcpy(alias, b, sizeof(alias));
			sub_mod_384x384(alias, a, alias, get_BLST_P().l);
			EXPECT_EQ(limbs_to_bytes(alias), expected_sub);
		}
	}
}

TEST(AssemblyArithmeticAudit, DivisionHelpersMatchBitwiseOracle)
{
	const limb_t divisor[2] = {
		TO_LIMB_T(0x0000000100000000),
		TO_LIMB_T(0xac45a4010001a402)
	};
	const Bytes<16> divisor_bytes = limbs_to_bytes(divisor);
	Bytes<16> one{};
	one[0] = 1;
	const std::array<Bytes<16>, 3> remainders = {
		Bytes<16>{}, one, decrement(divisor_bytes)
	};
	const limb_t low_words[] = {
		0, 1, TO_LIMB_T(0x0123456789abcdef), ~limb_t(0)
	};

	for (const Bytes<16>& remainder : remainders) {
		for (limb_t low : low_words) {
			limb_t dividend[3] = { low, 0, 0 };
			Bytes<24> dividend_bytes{};
			for (size_t i = 0; i < sizeof(low); ++i)
				dividend_bytes[i] = uint8_t(low >> (8 * i));
			for (size_t i = 0; i < remainder.size(); ++i)
				dividend_bytes[8 + i] = remainder[i];
			bytes_to_limbs(dividend, dividend_bytes);
			const auto expected = divide_bits(dividend_bytes, divisor_bytes);
			const limb_t exact_q = low_u64(expected.quotient);
			const limb_t estimate = div_3_limbs_v(dividend + 1,
				divisor[0], divisor[1]);
			EXPECT_TRUE(estimate == exact_q || estimate == exact_q + 1);

			limb_t corrected[3];
			std::memcpy(corrected, dividend, sizeof(corrected));
			EXPECT_EQ(quot_rem_128_v(corrected, divisor, estimate), exact_q);
			EXPECT_EQ(corrected[0], low_u64(expected.remainder));
			Bytes<16> expected_high = expected.remainder;
			for (size_t i = 0; i < 8; ++i)
				expected_high[i] = expected.remainder[8 + i];
			EXPECT_EQ(corrected[1], low_u64(expected_high));
			EXPECT_EQ(corrected[2], exact_q);
		}
	}

	EXPECT_EQ(div_3_limbs_v(divisor, divisor[0], divisor[1]), ~limb_t(0));

	const limb_t z[1] = { TO_LIMB_T(0xd201000000010000) };
	const Bytes<8> z_bytes = limbs_to_bytes(z);
	const limb_t high_words[] = { 0, 1, z[0] - 1 };
	for (limb_t high : high_words) {
		for (limb_t low : low_words) {
			limb_t dividend[2] = { low, high };
			const Bytes<16> dividend_bytes = limbs_to_bytes(dividend);
			const auto expected = divide_bits(dividend_bytes, z_bytes);
			const limb_t exact_q = low_u64(expected.quotient);
			limb_t corrected[2] = { low, high };
			EXPECT_EQ(quot_rem_64_v(corrected, z, exact_q), exact_q);
			EXPECT_EQ(corrected[0], low_u64(expected.remainder));
			EXPECT_EQ(corrected[1], exact_q);
		}
	}
}
#endif

TEST(AssemblyArithmeticAudit, SelectedZeroShiftsAreExactCopies)
{
	limb_t scalar[NLIMBS(256)] = { 1, 2, 3, 4 };
	limb_t scalar_out[NLIMBS(256)]{};
	lshift_mod_v256(scalar_out, scalar, 0, get_BLST_r().l);
	EXPECT_EQ(std::memcmp(scalar_out, scalar, sizeof(scalar)), 0);
	rshift_mod_v256(scalar_out, scalar, 0, get_BLST_r().l);
	EXPECT_EQ(std::memcmp(scalar_out, scalar, sizeof(scalar)), 0);

	limb_t field[NLIMBS(384)] = { 1, 2, 3, 4, 5, 6 };
	limb_t field_out[NLIMBS(384)]{};
	lshift_mod_v384(field_out, field, 0, get_BLST_P().l);
	EXPECT_EQ(std::memcmp(field_out, field, sizeof(field)), 0);
	rshift_mod_v384(field_out, field, 0, get_BLST_P().l);
	EXPECT_EQ(std::memcmp(field_out, field, sizeof(field)), 0);
}

#undef CATBLST_AUDIT_X86_64
