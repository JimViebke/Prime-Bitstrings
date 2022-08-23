#pragma once

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include <stdint.h>
#include <nmmintrin.h>

#include "config.hpp"
#include "boost_sqrt.hpp"
#include "util/types.hpp"

namespace mbp
{

	constexpr bool brute_force_is_prime(const size_t n)
	{
		if (n == 0) return false;
		if (n == 1) return false;

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
		using namespace mbp;

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

	template<size_t n_of_primes, size_t p_index = 1>
	__forceinline bool b2_has_small_divisor(const size_t number)
	{
		if constexpr (n_of_primes == p_index)
			return false;
		else
			return (number % small_primes_lookup[p_index] == 0) ? true : b2_has_small_divisor<n_of_primes, p_index + 1>(number);
	}

	bool mpir_is_prime(const mpz_class& p, gmp_randclass& r)
	{
		return bool(mpz_likely_prime_p(p.get_mpz_t(), r.get_randstate_t(), 0));
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
		/*
		The second arg to gcd() is a product of primes 3 through 13.
		This technically should be 2 through 13, but the alternating bitsum can't be 2 anyway:
		If a + b == (a prime number >11), then abs(a - b) is odd and never has 2 as a factor.
		*/

		size_t lookup = 0;

		// set bit i high if the GCD of (i, [product of primes]) is 1
		for (size_t i = 0; i < 32; ++i)
		{
			lookup |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << i;
		}

		return lookup;

		/*
		Alternatively, the loop can be:
			for (size_t i = 0; i < 32; ++i)
			{
				lookup |= size_t(gcd( ... ) == 1ull) << (32 - i);
				lookup |= size_t(gcd( ... ) == 1ull) << (32 + i);
			}

		Which changes the following:
			(lookup & (1ull << abs(pca - pcb))) == 0

		into this, requiring one less instruction:
			(lookup & (1ull << (pca + 32 - pcb))) == 0

		Any runtime difference was indistinguishable on my machine.
		*/
	}

	__forceinline void lex_permute(size_t& n)
	{
		const size_t t = n | (n - 1); // t gets v's least significant 0 bits set to 1
		// Next set to 1 the most significant bit to change, 
		// set to 0 the least significant ones, and add the necessary 1 bits.

		unsigned long idx;
		_BitScanForward64(&idx, n);

		// suppress "warning C4146 : unary minus operator applied to unsigned type, result still unsigned"
	#pragma warning (push)
	#pragma warning (disable: 4146)
		n = (t + 1) | (((~t & -~t) - 1) >> (idx + 1));
	#pragma warning (pop)
	}

}
