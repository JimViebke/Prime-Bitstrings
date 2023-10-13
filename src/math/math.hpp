#pragma once

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include <bit>
#include <immintrin.h>
#include <stdint.h>
#include <vector>

#include "../config.hpp"
#include "../util/types.hpp"

namespace mbp
{
	namespace detail
	{
		constexpr std::vector<sieve_prime_t> build_small_primes_lookup_impl()
		{
			std::vector<sieve_prime_t> primes;
			primes.push_back(2);

			// capture requirements for the small primes lookup here
			constexpr size_t largest_prime = std::max({ prime_sieve::largest_sieve_prime,
													  prime_test::mpir_trial_div_cap });
			std::vector<uint8_t> sieve(largest_prime + 1, true);

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

	inline bool mpir_is_prime(const mpz_class& p, gmp_randclass& r)
	{
		return mpz_likely_prime_p(p.get_mpz_t(), r.get_randstate_t(), 0);
	}

	constexpr __forceinline auto pop_count(uint64_t n)
	{
		if (std::is_constant_evaluated())
		{
			using T = decltype(_mm_popcnt_u64(n));
			return (T)std::popcount(n);
		}
		else
		{
			return _mm_popcnt_u64(n);
		}
	}

	consteval uint64_t build_tiny_primes_lookup()
	{
		// Generate a 64 bit lookup, where the prime-numbered bits are set high
		// 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61

		uint64_t lookup = 0;

		// The popcount of a p13 can't be divisible by 2, 3, 5, 7 or 11
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

}
