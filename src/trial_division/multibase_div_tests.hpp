#pragma once

/*
Determine if a bitstring has a small prime divisor in base b.

Given any bitstring, the value of each digit is either base^(digit position) or 0. For example, in base 7, bits 0, 1, and 2 can only represent 7^0, 7^1, 7^2... or 0.

Therefore, for base b, place value n, and a small prime p, we calculate and store b^n % p. These remainders eventually follow a repeating pattern, often shorter than the chosen n, so we only need to store k unique remainders. At runtime, for each of the k remainders, use a bitmask to select every k-th bit, perform a popcount, and multiply the count by that k-th remainder. Adding these k results together gives us a small sum of remainders, s. If s is evenly divisible by p, then the bitstring is divisible by p in base b. Discard it.
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

					// Skip always-composite cases
					if (base % p == 0) continue;

					// Always suppress hardcoded div tests

					// Hardcoded div tests with 4 remainders
					if (base == 3 && p == 5) continue;
					if (base == 5 && p == 13) continue;
					if (base == 8 && p == 13) continue;
					if (base == 4 && p == 17) continue;

					// 3 remainders
					if (base == 4 && p == 7) continue;
					if (base == 3 && p == 13) continue;
					if (base == 9 && p == 13) continue;

					// 6 remainders
					if (base == 3 && p == 7) continue;
					if (base == 5 && p == 7) continue;
					if (base == 4 && p == 13) continue;
					if (base == 10 && p == 13) continue;

					// 5 remainders
					if (base == 3 && p == 11) continue;
					if (base == 4 && p == 11) continue;
					if (base == 5 && p == 11) continue;
					if (base == 9 && p == 11) continue;

					// 10 remainders
					if (base == 6 && p == 11) continue;
					if (base == 7 && p == 11) continue;
					if (base == 8 && p == 11) continue;

					// 12 remainders
					if (base == 6 && p == 13) continue;
					if (base == 7 && p == 13) continue;
					if (base == 11 && p == 13) continue;

					// 16 remainders
					if (base == 3 && p == 17) continue;
					if (base == 5 && p == 17) continue;
					if (base == 6 && p == 17) continue;
					if (base == 7 && p == 17) continue;
					if (base == 10 && p == 17) continue;
					if (base == 11 && p == 17) continue;
					if (base == 12 && p == 17) continue;

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

				#endif

					uncompressed_div_test_t dt{ .base = base_t(base), .prime_idx = prime_idx_t(i) };

					// calculate base^j mod prime, where j is the place value
					for (size_t j = 0; j < 64; ++j)
					{
						remainder_t rem = remainder_t(pk::powMod(base, j, p));
						if (rem == 1 && j > 0)
						{
							// The pattern is repeating; stop generating further terms
							break;
						}

						dt.remainders[j] = rem;
						dt.n_of_remainders++;
					}

					uncompressed_dts.push_back(dt);
				}
			}

			// Div tests are periodically re-ordered during runtime based on performance.
			// To reflect real-world performance in benchmarks, pre-order the div tests based
			// on cached performance data (the hitcount of each test).
			// Otherwise, start with the default ordering of primes, low to high.
			if constexpr (benchmark_mode)
			{
				for (auto& dt : uncompressed_dts)
				{
					dt.hits = cached_hitcount_for(dt.base, small_primes_lookup[dt.prime_idx]);
				}

				std::sort(uncompressed_dts.begin(), uncompressed_dts.end(),
						  [](const auto& a, const auto& b) { return a.hits > b.hits; });
			}

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

		// calculate a^b mod m
		template<size_t a, size_t b, size_t m>
		struct pow_mod
		{
			static constexpr size_t rem = pk::powMod(a, b, m);
		};

		// bitmask for all base^n mod prime
		template<size_t base, size_t prime>
		struct bitmask_for
		{
			static consteval size_t f()
			{
				// Calculate the period of base^n % prime
				size_t i = 0;
				for (; i < 64; ++i)
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
	static div_tests_t div_tests = generate_div_tests(); // intellisense false positive



	namespace detail
	{
		consteval size_t calculate_prime_factor_lookup_size()
		{
			constexpr div_tests_t div_tests_constexpr = generate_div_tests();

			size_t largest_sum = 0;

			// Calculate the largest possible sum of remainders of a number with every bit set
			for (const auto& div_test : div_tests_constexpr)
			{
				const size_t sum = std::accumulate(div_test.remainders.begin(), div_test.remainders.end(), size_t(0));

				if (sum > largest_sum)
					largest_sum = sum;
			}

			// + 1 so "lookup[sum]" is always in range
			return largest_sum + 1;
		}

		using prime_lookup_t = util::narrowest_uint_for_n_bits<div_test::n_of_primes>;
		constexpr size_t prime_factor_lookup_size = calculate_prime_factor_lookup_size(); // intellisense false positive

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
					lookup[j] |= (0b1 << i);
				}
			}

			return lookup;
		}

		const std::vector<prime_lookup_t> prime_factor_lookup = build_prime_factor_lookup_old();
	}

	// Replaces "n % prime[idx] == 0" with "lookup[n] & (1 << idx)", usually as bittest + cmov
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
