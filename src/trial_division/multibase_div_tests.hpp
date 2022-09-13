#pragma once

/*
To determine if a bitstring has a small prime divisor in base b, we can do better than converting the bitstring to base b, then calculating bistring % smallprime == 0.

Given any bitstring in the original base 2, the value of each digit is either 0 or base^(digit position). For example, in base 7, the bits can only represent 7^0, 7^1, 7^2...

Therefore, for base b, place value n, and a small prime p, we can easily calculate and store b^n % p. These remainders follow a repeating pattern often shorter than the chosen n, so we only need to store k unique remainders. At runtime, for each of the k remainders, use a bitmask to select every k-th bit, run a popcount, and multiply by the k-th remainder. Adding these k results together gives us a small sum of remainders, s. If s is evenly divisible by p, the bitstring is divisible by p in base b. Discard it.
*/

#include <vector>

#include "../math/math.hpp"
#include "../math/pk_prime.hpp"
#include "types.hpp"
#include "div_test_hits.hpp"

namespace mbp::div_test
{
	namespace detail
	{
		consteval std::vector<div_test_t> generate_div_tests_impl()
		{
			std::vector<uncompressed_div_test_t> uncompressed_dts;

			for (size_t i = 1; i < n_of_primes; ++i) // for each small prime starting from 3
			{
				for (size_t base = 3; base <= up_to_base; ++base) // for each base 3..n
				{
					const auto p = small_primes_lookup[i];

					// Always suppress hardcoded div tests
					if (base == 3 && p == 5) continue;

					if (base == 4 && p == 7) continue;

					// Hardcoded div tests with 6 remainders
					if (base == 3 && p == 7) continue;
					if (base == 5 && p == 7) continue;
					if (base == 4 && p == 13) continue;
					if (base == 10 && p == 13) continue;

					// another block of hardcoded tests, ordered by their
					// original measured hit counts
					if (base == 8 && p == 13) continue;
					if (base == 5 && p == 13) continue;
					if (base == 4 && p == 17) continue;

					// Hardcoded 5-remainder div tests by 11
					if (base == 3 && p == 11) continue;
					if (base == 4 && p == 11) continue;
					if (base == 5 && p == 11) continue;
					if (base == 9 && p == 11) continue;

					// Hardcoded 10-remainder div tests by 11
					if (base == 6 && p == 11) continue;
					if (base == 7 && p == 11) continue;
					if (base == 8 && p == 11) continue;

				#if !analyze_div_tests or suppress_extra_div_tests
					if (base == 4 && p == 3) continue; //  base  4^n % 3 unused
					if (base == 5 && p == 3) continue; //  base  5^n % 3 unused
					if (base == 7 && p == 3) continue; //  base  7^n % 3 unused
					if (base == 8 && p == 3) continue; //  base  8^n % 3 unused
					if (base == 10 && p == 3) continue; // base 10^n % 3 unused
					if (base == 11 && p == 3) continue; // base 11^n % 3 unused

					if (base == 4 && p == 5) continue; //  base  4^n % 5 unused
					if (base == 6 && p == 5) continue; //  base  6^n % 5 unused
					if (base == 7 && p == 5) continue; //  base  7^n % 5 unused
					if (base == 9 && p == 5) continue; //  base  9^n % 5 unused
					if (base == 11 && p == 5) continue; // base 11^n % 5 unused
					if (base == 12 && p == 5) continue; // base 12^n % 5 unused

					if (base == 6 && p == 7) continue; //  base  6^n % 7 unused
					if (base == 8 && p == 7) continue; //  base  8^n % 7 unused
					if (base == 9 && p == 7) continue; //  base  9^n % 7 unused

					if (base == 10 && p == 11) continue; // base 10^n % 11 unused
					if (base == 12 && p == 11) continue; // base 12^n % 11 unused

					if (base == 12 && p == 13) continue; // base 12^n % 13 unused

					// If two div tests are effectively identical, remove one
					if (base == 8 && p == 5) continue; //  base  8^n % 5 is congruent to 3^n % 5
					if (base == 10 && p == 7) continue; // base 10^n % 7 is congruent to 3^n % 7
					if (base == 11 && p == 7) continue; // base 11^n % 7 is congruent to 4^n % 7 
					if (base == 12 && p == 7) continue; // base 12^n % 7 is congruent to 5^n % 7

					// removed for being unused (due to ordering)
					if (base == 9 && p == 73) continue; //   base  9^n %  73:   -       6 remainders : 1   9   8  72  64  65
				#endif

					uncompressed_div_test_t dt{ .base = base_t(base), .prime_idx = prime_idx_t(i) };

					// calculate base^j mod prime, where j is the place value
					for (size_t j = 0; j < max_remainders; ++j)
					{
						remainder_t rem = remainder_t(pk::powMod(base, j, small_primes_lookup[i]));
						if (rem == 1 && j > 0)
						{
							// The pattern is repeating - store what we have, then break
							uncompressed_dts.push_back(dt);
							break;
						}

						dt.remainders[j] = rem;
						dt.n_of_remainders++;
					}

					// Special case where we maxed out our terms without finding a repeat.
					// Save this div test if and only if the next term is 1 (ie, a repeat)
					if (dt.n_of_remainders == max_remainders &&
						pk::powMod(base, max_remainders, small_primes_lookup[i]) == 1)
					{
						uncompressed_dts.push_back(dt);
					}
				}
			}

			for (auto& dt : uncompressed_dts)
			{
				dt.hits = cached_hitcount_for(dt.base, small_primes_lookup[dt.prime_idx]);
			}

			// Order div tests by prime divisor
			std::sort(uncompressed_dts.begin(), uncompressed_dts.end(), [](const auto& a, const auto& b)
					  {
						  //if (a.prime_idx == b.prime_idx)
							 // return a.n_of_remainders < b.n_of_remainders;
						  // return a.prime_idx < b.prime_idx;

						  //if (a.n_of_remainders == b.n_of_remainders)
							 // return a.prime_idx < b.prime_idx;
						  //return a.n_of_remainders < b.n_of_remainders;

						  //return
							 // size_t(a.n_of_remainders) * small_primes_lookup[a.prime_idx] <
							 // size_t(b.n_of_remainders) * small_primes_lookup[b.prime_idx];

						  return
							  1.0 / double(a.hits) <
							  1.0 / double(b.hits);

						  //return
							 // double(a.n_of_remainders) * double(small_primes_lookup[a.prime_idx]) / (1. * double(a.hits)) <
							 // double(b.n_of_remainders) * double(small_primes_lookup[b.prime_idx]) / (1. * double(b.hits));
					  });

			std::vector<div_test_t> div_tests;
			div_tests.reserve(uncompressed_dts.size());

			for (const auto& udt : uncompressed_dts)
			{
				div_test_t dt{
					.prime_idx = udt.prime_idx,
					.n_of_remainders = udt.n_of_remainders,
					.remainders = udt.remainders };

				// Repeat the remainders so all div tests have 64 terms
				for (size_t i = dt.n_of_remainders; i < 64; ++i)
				{
					dt.remainders[i] = dt.remainders[i - dt.n_of_remainders];
				}

			#if analyze_div_tests
				// We do need to copy this if we're in analyze mode
				dt.base = udt.base;
			#endif

				div_tests.push_back(dt);
			}

			return div_tests;
		}

		// bitmask for all base^n mod prime
		template<size_t base, size_t prime>
		struct bitmask_for
		{
			static consteval size_t f()
			{
				size_t i = 0;
				for (; i < 64; ++i) // calculate base^i MOD prime
				{
					remainder_t rem = remainder_t(pk::powMod(base, i, prime));
					if (rem == 1 && i > 0) break; // break when the pattern repeats
				}

				size_t bitmask = 0;
				for (size_t j = 0; j < 64 && i < 64; j += i)
				{
					bitmask <<= i;
					bitmask |= 1;
				}

				return bitmask;
			}
			static constexpr size_t val = f();
		};

		// calculate a^b mod m
		template<size_t a, size_t b, size_t m>
		struct pow_mod
		{
			static constexpr size_t rem = pk::powMod(a, b, m);
		};

		template<size_t bitmask>
		struct period_of
		{
			static consteval size_t f()
			{
				for (size_t i = 1; i < 64; ++i)
					if ((bitmask >> i) & 1)
						return i;
			}
			static constexpr size_t val = f();
		};

		template<size_t prime>
		struct get_prime_index
		{
			static consteval size_t f()
			{
				for (size_t i = 0; i < 64; ++i)
					if (small_primes_lookup[i] == prime)
						return i;
			}
			static constexpr size_t idx = f();
		};
	}



	// looping, sorted div tests:

	constexpr size_t div_tests_size = detail::generate_div_tests_impl().size();
	consteval std::array<div_test_t, div_tests_size> generate_div_tests()
	{
		std::array<div_test_t, div_tests_size> div_tests{};
		const auto x = detail::generate_div_tests_impl();
		std::copy(x.begin(), x.end(), div_tests.begin());
		return div_tests;
	}

	using div_tests_t = std::array<div_test::div_test_t, div_test::div_tests_size>;
	static div_test_constexpr div_tests_t div_tests = generate_div_tests(); // intellisense false positive



	namespace detail
	{
		consteval size_t calculate_prime_factor_lookup_size()
		{
			size_t largest_remainders_sum = 0;

		#if analyze_div_tests
			// Gross, but if we're in "analyze" mode, we'll need to generate our own constexpr version of div_tests
		#pragma warning (push)
		#pragma warning (disable: 4459)
			constexpr div_tests_t div_tests = generate_div_tests();
		#pragma warning (pop)
		#endif

			// for each div test
			for (const auto& div_test : div_tests)
			{
				// calculate the sum of having every bit set
				size_t remainders_sum = 0;
				for (size_t i = 0; i < 64; ++i)
					remainders_sum += div_test.remainders[i % div_test.n_of_remainders];

				// keep track of the largest sum
				if (remainders_sum > largest_remainders_sum)
					largest_remainders_sum = remainders_sum;
			}

			// + 1 so "lookup[sum]" is always safe.
			return largest_remainders_sum + 1;
		}

		using prime_lookup_t = narrowest_uint_for_n_bits<div_test::n_of_primes>;
		constexpr size_t prime_factor_lookup_size = calculate_prime_factor_lookup_size();

		// Faster version
		std::vector<prime_lookup_t> build_prime_factor_lookup_old()
		{
			std::vector<prime_lookup_t> lookup;
			lookup.reserve(prime_factor_lookup_size);

			for (size_t i = 0; i < prime_factor_lookup_size; ++i)
			{
				prime_lookup_t entry = 0;
				for (prime_lookup_t p = 0; p < div_test::n_of_primes; ++p)
				{
					entry |= (prime_lookup_t((i % small_primes_lookup[p]) == 0) << p);
				}

				lookup.push_back(entry);
			}

			return lookup;
		}

		// Slower version
		std::vector<prime_lookup_t> build_prime_factor_lookup_new()
		{
			std::vector<prime_lookup_t> lookup(prime_factor_lookup_size, 0);

			for (size_t i = 0; i < div_test::n_of_primes; ++i)
			{
				const prime_lookup_t p = small_primes_lookup[i];

				for (prime_lookup_t j = 0; j < prime_factor_lookup_size; j += p)
				{
					lookup[j] |= (0x1 << i);
				}
			}

			return lookup;
		}

		// Replaces "n % prime[k] == 0" with "lookup[n] & (1 << k)"
		const std::vector<prime_lookup_t> prime_factor_lookup = build_prime_factor_lookup_old();
	}

	__forceinline bool has_small_prime_factor(const size_t n, const prime_idx_t prime_index)
	{
		return (detail::prime_factor_lookup[n] & (detail::prime_lookup_t(1) << prime_index)) != 0;
	}



	// unrolled div tests using recursive templates:

	template<size_t base>
	struct in_base
	{
		static constexpr size_t val = base;
	};

	// Suppress warnings about bitmasks having upper bits moved
#pragma warning(push)
#pragma warning(disable: 26450)

	namespace detail
	{
		template<size_t divisor, size_t base, size_t mask, size_t place_value>
		__forceinline void recursive_is_divisible_by(size_t& rem, const size_t number)
		{
			rem += pop_count(number & (mask << place_value)) * pow_mod<base, place_value, divisor>::rem;
			if constexpr (place_value > 0)
				recursive_is_divisible_by<divisor, base, mask, place_value - 1>(rem, number);
		}
	}

	template<size_t divisor, typename base_t>
	__forceinline bool recursive_is_divisible_by(const size_t number)
	{
		using namespace detail;

		constexpr size_t base = base_t::val;
		static_assert(base >= 3 && base <= up_to_base);
		constexpr size_t bitmask = bitmask_for<base, divisor>::val;
		size_t rem = 0;

		recursive_is_divisible_by<divisor, base, bitmask, period_of<bitmask>::val - 1>(rem, number);

		return has_small_prime_factor(rem, get_prime_index<divisor>::idx);
	}

#pragma warning(pop)
}
