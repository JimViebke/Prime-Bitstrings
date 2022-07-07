#pragma once

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include <stdint.h>
#include <nmmintrin.h>

#include "config.hpp"
#include "franken_boost.hpp"

constexpr bool brute_force_is_prime(const size_t n)
{
	if (n == 2) return true;

	if (n % 2 == 0) return false;

	const size_t sqrt_n = size_t(franken_boost::sqrt(n));
	for (size_t i = 3; i <= sqrt_n; i += 2)
	{
		if (n % i == 0) return false;
	}

	return true;
}

namespace detail
{
	constexpr const std::vector<size_t> build_small_primes_lookup_impl()
	{
		std::vector<size_t> primes;

		for (size_t i = 2; i <= mbp::sieve_primes_cap; ++i)
			if (brute_force_is_prime(i))
				primes.push_back(i);

		return primes;
	}

	constexpr std::array<size_t, build_small_primes_lookup_impl().size()> build_small_primes_lookup()
	{
		decltype(build_small_primes_lookup()) primes{};
		const auto x = build_small_primes_lookup_impl();
		std::copy(x.begin(), x.end(), primes.begin());
		return primes;
	}
}

constexpr std::array small_primes_lookup = detail::build_small_primes_lookup();

namespace gmp_random
{
	gmp_randclass r{ gmp_randinit_mt };
}

bool mpir_is_prime(const mpz_class& p, const size_t div = 0)
{
	return bool(mpz_likely_prime_p(p.get_mpz_t(), gmp_random::r.get_randstate_t(), div));
}

size_t binary_to_base(size_t binary, const size_t base)
{
	constexpr size_t highest_bit = size_t(1) << 63;

	size_t num = 0;

	// for each bit
	for (auto i = 0; i < 64; binary <<= 1, ++i)
	{
		num *= base;

		// if the (i-1)th bit is set, +1
		if ((binary & highest_bit) == highest_bit)
		{
			num += 1;
		}
	}

	return num;
}

mpz_class bin_to_base(const mpz_class& binary, const int base)
{
	return mpz_class{ binary.get_str(2), base };
}

inline auto pop_count(uint64_t n)
{
	return _mm_popcnt_u64(n);
}

constexpr uint64_t build_tiny_primes_lookup()
{
	// Generate a 64 bit lookup, where the prime-numbered bits are set high
	// 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61

	uint64_t lookup = 0;

	// We don't care about tiny pNs, therefore leave out the smallest of the primes.
	//lookup |= (1ull << 2);
	//lookup |= (1ull << 3);
	//lookup |= (1ull << 5);
	//lookup |= (1ull << 7);
	//lookup |= (1ull << 11);
	lookup |= (1ull << 13);
	lookup |= (1ull << 17);
	lookup |= (1ull << 19);
	lookup |= (1ull << 23);
	lookup |= (1ull << 29);
	lookup |= (1ull << 31);
	lookup |= (1ull << 37);
	lookup |= (1ull << 41);
	lookup |= (1ull << 43);
	lookup |= (1ull << 47);
	lookup |= (1ull << 53);
	lookup |= (1ull << 59);
	lookup |= (1ull << 61);

	return lookup;
}

constexpr size_t gcd(size_t a, size_t b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);
}

constexpr size_t build_gcd_1155_lookup()
{
	size_t lookup = 0;

	// set bit i high if and only if the GCD of (i, 1155) is 1
	for (size_t i = 0; i < 32; ++i)
	{
		lookup |= size_t(gcd(i, 1155u) == 1ull) << i;
	}

	return lookup;
}
