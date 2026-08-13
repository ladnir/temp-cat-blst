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
#include "blst/fields/Fr.h"

using namespace blst;

namespace
{
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
