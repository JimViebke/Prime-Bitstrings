#pragma once

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

// Hacky way to hide the contents of zmmintrin.h
// MSVC always includes it from immintrin.h
#ifndef _ZMMINTRIN_H_INCLUDED
#define _ZMMINTRIN_H_INCLUDED
#endif

#include <immintrin.h>
#include <stdint.h>

#include "config.hpp"
#include "util/types.hpp"

namespace mbp
{
	namespace detail
	{
		consteval std::vector<sieve_prime_t> build_small_primes_lookup_impl()
		{
			std::vector<sieve_prime_t> primes;
			primes.push_back(2);

			std::vector<uint8_t> sieve(sieve_primes_cap + 1, true);

			for (size_t i = 3; i < sieve.size(); i += 2)
			{
				if (!sieve[i]) continue;

				primes.push_back(sieve_prime_t(i));

				// mark off all odd multiples of i, except for i
				for (size_t j = 3 * i; j < sieve.size(); j += 2 * i)
				{
					sieve[j] = false;
				}
			}

			return primes;
		}

		constexpr size_t n_of_small_primes = build_small_primes_lookup_impl().size();
		consteval std::array<sieve_prime_t, n_of_small_primes> build_small_primes_lookup()
		{
			decltype(build_small_primes_lookup()) primes{};
			const auto x = build_small_primes_lookup_impl();
			std::copy(x.begin(), x.end(), primes.begin());
			return primes;
		}
	}

	constexpr std::array small_primes_lookup = detail::build_small_primes_lookup();

	bool mpir_is_prime(const mpz_class& p, gmp_randclass& r)
	{
		return mpz_likely_prime_p(p.get_mpz_t(), r.get_randstate_t(), 0);
	}

	mpz_class bin_to_base(const mpz_class& binary, const size_t base)
	{
		return mpz_class{ binary.get_str(2), int(base) };
	}

	inline auto pop_count(uint64_t n)
	{
		return _mm_popcnt_u64(n);
	}

	consteval uint64_t build_tiny_primes_lookup()
	{
		// Generate a 64 bit lookup, where the prime-numbered bits are set high
		// 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61

		uint64_t lookup = 0;

		// The popcount of a p12 can't be divisible by 2, 3, 5, 7 or 11
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

	consteval size_t build_gcd_lookup()
	{
		// The second arg to gcd() is a product of primes 3 through 13.
		// This should be 2 through 13, but the alternating bitsum can't be 2:
		// If a + b == (a prime number >11), then abs(a - b) is odd and never has 2 as a factor.

		size_t lookup = 0;

		// set bit i high if the GCD of (i, [product of primes]) is 1
		for (size_t i = 0; i < 32; ++i)
		{
			lookup |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << (32 - i);
			lookup |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << (32 + i);
		}

		return lookup;

		/*
		Loop impl 1:
			lookup |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << i;

		Use:
			(lookup & (1ull << abs(pca - pcb))) == 0

		Loop impl 2:
			lookup |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << (32 - i);
			lookup |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << (32 + i);

		Use:
			(lookup & (1ull << (pca + 32 - pcb))) == 0
		*/
	}

}
