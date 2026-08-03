
#pragma once

#include "blst/common/defines.h"
#include <array>
#include <vector>
#include <random>
#include <span>
#include "blst/fields/Fr.h"

namespace blst
{


	struct Factor {
		Factor() = default;
		Factor(const Factor&) = default;
		Factor(u64 f, u64 e) :mFactor(f), mExp(e) {}

		u64 mFactor;
		u64 mExp;
	};

	inline std::vector<Factor> uniqueFactor(u64 x)
	{
		std::vector<Factor> r;
		u64 X = x;
		if (X < 2) return r;

		// Factor out powers of 2
		if ((X & 1) == 0)
		{
			u64 e = 0;
			do { X >>= 1; ++e; } while ((X & 1) == 0);
			r.emplace_back(2, e);
		}

		// Trial divide by odd numbers up to sqrt(X)
		for (u64 i = 3; i <= X / i; i += 2)
		{
			if (X % i == 0)
			{
				u64 e = 0;
				do { X /= i; ++e; } while (X % i == 0);
				r.emplace_back(i, e);
			}
		}

		// Any remaining factor is prime
		if (X > 1)
			r.emplace_back(X, 1);

		return r;
	}

	template<typename F>
	F findGenerator() {
		auto p = F::order();

		// we are interested in the range [2, p-2]
		auto dist = [&]() {
			F f;
			do {
				f = F::rand();
			} while (f < 2 || f >= p - 2);

			return f;
			};
		auto factors = uniqueFactor(p - 1);

		while (true) {
			F a = dist();
			bool ok = true;
			for (auto q : factors)
			{
				if (a.pow((p - 1) / q.mFactor) == 1)
				{
					ok = false;
					break;
				}
			}

			if (ok)
				return a;
		}
	}



	// returns true if u is a primitive root of unity.
	// factors should be the unique factors of n.
	template<typename F>
	inline bool isPrimRootOfUnity(std::span<Factor> factors, const F& u)
	{
		auto p = F::order();
		u64 n = 1;
		for (auto fe : factors)
		{
			for (u64 i = 0; i < fe.mExp; ++i)
				n *= fe.mFactor;
		}

		if ((p - 1) % n)
			return false;

		// make suer u is in BasicFp*
		if (u == 0 || u.integer() % p == 0)
			return false;

		// make sure u is a root of unity.
		if (u.pow(n) != 1)
			return false;

		// check that u is a primitive root of unity.
		for (u64 i = 0; i < factors.size(); ++i)
		{
			if (u.pow(n / factors[i].mFactor) == 1)
			{
				return false;
			}
		}
		return true;
	}
	// returns true if f is a primitive root of unity.
	template<typename F>
	inline bool isPrimRootOfUnity(u64 n, const F& u)
	{
		auto factors = uniqueFactor(n);
		return isPrimRootOfUnity<F>(factors, u);
	}




	template<>
	inline bool isPrimRootOfUnity<Fr>(std::span<Factor> factors, const Fr& u)
	{
		u64 n = 1;
		for (auto fe : factors)
		{
			for (u64 i = 0; i < fe.mExp; ++i)
				n *= fe.mFactor;
		}

		// make suer u is in Fp*
		if (u == Fr::zero())
			return false;

		Fr uPowN;
		Fr::exp(uPowN, u, (byte*)&n, std::bit_width(n));
		// make sure u is a root of unity.
		if (uPowN != Fr::one())
			return false;

		// check that u is a primitive root of unity.
		for (u64 i = 0; i < factors.size(); ++i)
		{
			auto ni = n / factors[i].mFactor;
			Fr uPowNi;
			Fr::exp(uPowNi, u, (byte*)&ni, std::bit_width(ni));
			if (uPowNi == Fr::one())
			{
				return false;
			}
		}
		return true;
	}

	// returns true if f is a primitive root of unity.
	template<>
	inline bool isPrimRootOfUnity<Fr>(u64 n, const Fr& u)
	{
		auto factors = uniqueFactor(n);
		return isPrimRootOfUnity<Fr>(factors, u);
	}


	// return a primitive n-root of unity given a generator.
	template<typename F>
	inline F primRootOfUnity(u64 n, F generator)
	{
		auto p = F::order();
		return generator.pow((p - 1) / n);
	}

	template<typename F>
	inline F primRootOfUnity(u64 n)
	{
		F G = findGenerator<F>();
		return primRootOfUnity<F>(n, G);
	}



	// Explicit specialization of primRootOfUnity for Fr.
	// This bypasses the random search in findGenerator by using the known
	// structure of the BLS12-381 scalar field (2-adicity = 32, generator = 7).
	template<>
	inline Fr primRootOfUnity<Fr>(u64 n)
	{
		// The explicit table cache is initialized exactly once (thread-safe).
		static const std::array<Fr, 33> roots = []() {
			std::array<Fr, 33> table; // Stores primitive 2^k-th roots for k=0..32

			table[32][0] = 0xb9b58d8c5f0e466a;
			table[32][1] = 0x5b1b4c801819d7ec;
			table[32][2] = 0x0af53ae352a31e64;
			table[32][3] = 0x5bf3adda19e9b27b;

			// Compute 2^k roots by repeated squaring
			for (int i = 31; i >= 0; --i) {
				Fr::sqr(table[i], table[i + 1]);
			}
			return table;

			}();

		// Verify n is a power of 2
		u64 k = 0;
		if (n > 0 && (n & (n - 1)) == 0) {
			u64 v = n;
			while (v >>= 1) k++;
		}
		else {
			throw std::runtime_error("primRootOfUnity: n must be a power of 2");
		}

		if (k >= roots.size()) {
			throw std::runtime_error("primRootOfUnity: n exceeds 2-adicity (2^32)");
		}

		// roots[k] is the primitive 2^k-th root of unity
		return roots[k];
	}


	// returns true if u is an n-root of unity.
	template<typename F>
	inline bool isRootOfUnity(u64 n, const F& u)
	{
		auto p = F::order();
		if (u == 0)
			return false;
		if (u.pow(n) != 1)
			return false;
		return true;
	}
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define XX_RTE_LOC std::runtime_error(__FILE__ ":" STR(__LINE__))


	// Pre-compute powers of a primitive n-th root of unity.
	template<typename F>
	void precomputeRootsOfUnity(const F& psi, std::span<F> roots)
	{
		auto n = roots.size();
		if (n == 0)
			throw XX_RTE_LOC;
		auto ln = log2ceil(n);
		auto qq = F::order() - 1;
		auto q = static_cast<std::conditional_t<(sizeof(decltype(qq)) < sizeof(u64)), u64, decltype(qq)>>(qq);
		if (n != 1ull << ln)
			throw XX_RTE_LOC;
		if (n > q)
			throw XX_RTE_LOC;
		if (q % n != 0)
			throw XX_RTE_LOC;
		if (psi.pow(n) != 1)
			throw XX_RTE_LOC;
#ifndef NDEBUG
		if (isPrimRootOfUnity(n, psi) == false)
			throw XX_RTE_LOC;
#endif

		F root = F(1);
		for (u64 i = 0; i < n; ++i)
		{
			roots[i] = root;
			root *= psi;
		}
	}



}
