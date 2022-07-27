#pragma once

#include <vector>

#include "math/math.hpp"
#include "math/pk_prime.hpp"

/*
To determine if a bitstring has a small prime divisor in base b, we can do better than converting the bitstring to base b, then calculating bistring % smallprime == 0.

Given any bitstring in the original base 2, the value of each digit is either 0 or base^(digit position). For example, in base 3, the bits can only represent 3^0, 3^1, 3^2... respectively.

Therefore, we can easily calculate the remainder of [place value] % [smallprime]. This could result in a 64-entry lookup, however, these remainders follow a short, repetitive pattern. Therefore, instead of comparing 64 bits with 64 different remainders, we only need to grab each unique remainder. Given k unique remainders, we can use a mask k times to select every k-th bit in the original number, perform a popcount, and multiple that popcount by the remainder. Adding these k remainders together gives us a small number, not greater than perhaps 100-200.

Then, cheaply determine if this sum is divisible by p. This determines if the bitstring would be divisible by p in base b.
*/

namespace mbp::div_test
{
	using base_t = narrowest_uint_for_val<up_to_base>;
	using prime_idx_t = narrowest_uint_for_val<div_test::n_of_primes>;
	using n_of_remainders_t = narrowest_uint_for_val<div_test::max_remainders>;
	using remainder_t = sieve_prime_t;

	class div_test_t
	{
	public:
#if analyze_div_tests
		bool used = false;
		base_t base = 0;
		uint32_t hits = 0;
#endif
		prime_idx_t prime_idx = 0;
		n_of_remainders_t n_of_remainders = 0; // is also the index of the req'd bitmask

#if USE_UNCACHED
#else
		bool is_first_with_n_remainders = false;
#endif

		std::array<remainder_t, max_remainders> remainders{ 0 };
	};

	namespace detail
	{
		consteval std::vector<div_test_t> generate_div_tests_impl()
		{
			std::vector<div_test_t> div_tests;

			for (size_t i = 1; i < n_of_primes; ++i) // for each small prime starting from 3
			{
				for (size_t base = 3; base <= up_to_base; ++base) // for each base 3..n
				{
					const auto p = small_primes_lookup[i];

					// Always suppress hardcoded div tests
					if (base == 3 && p == 5) continue;

					if (base == 3 && p == 7) continue;
					if (base == 4 && p == 7) continue;
					if (base == 5 && p == 7) continue;

#if !analyze_div_tests or suppress_extra_div_tests
					if (base == 2 && p == 3) continue;

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

					// removed for performance (may be due to ordering)
					// if (base == 8 && p == 73) continue; //   base  8^n %  73:  11 hits  3 remainders : 1   8  64
					// if (base == 7 && p == 43) continue; //   base  7^n %  43: 155 hits  6 remainders : 1   7   6  42  36  37
					// if (base == 10 && p == 101) continue; // base 10^n % 101:   4 hits  4 remainders : 1  10 100  91

					// removed for being unused (due to ordering)
					if (base == 9 && p == 73) continue; //   base  9^n %  73:   -       6 remainders : 1   9   8  72  64  65
#endif

#if analyze_div_tests
					div_test_t dt{ .base = base_t(base), .prime_idx = prime_idx_t(i) };
#else
					div_test_t dt{ .prime_idx = prime_idx_t(i) };
#endif

					// calculate base^j mod prime, where j is the place value
					for (size_t j = 0; j < max_remainders; ++j)
					{
						remainder_t rem = remainder_t(pk::powMod(base, j, small_primes_lookup[i]));
						if (rem == 1 && j > 0)
						{
							// The pattern is repeating - store what we have, then break
							div_tests.push_back(dt);
							break;
						}

						dt.remainders[j] = rem;
						dt.n_of_remainders++;
					}
				}
			}

			// Order div tests by worthwhileness
			std::sort(div_tests.begin(), div_tests.end(), [] (const auto& a, const auto& b)
					  {
						  //if (a.prime_idx == b.prime_idx)
							 // return a.n_of_remainders < b.n_of_remainders;
						  //return a.prime_idx < b.prime_idx;

						  //if (a.n_of_remainders == b.n_of_remainders)
							 // return a.prime_idx < b.prime_idx;
						  //return a.n_of_remainders < b.n_of_remainders;

						  return
							  size_t(a.n_of_remainders) * small_primes_lookup[a.prime_idx] <
							  size_t(b.n_of_remainders) * small_primes_lookup[b.prime_idx];
					  });

#if USE_UNCACHED
#else
			// Mark each div test that is the first test with N remainders
			for (size_t i = 0; i <= max_remainders; ++i)
			{
				for (auto& div_test : div_tests)
				{
					if (div_test.n_of_remainders == i)
					{
						div_test.is_first_with_n_remainders = true;
						break;
					}
				}
			}
#endif

			//for (auto& dt : div_tests)
			//	std::cout << "b" << size_t(dt.base) << " % " << small_primes_lookup[dt.prime_idx] << " \t" << size_t(dt.n_of_terms) << " terms\n";

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

		consteval size_t calculate_prime_factor_lookup_size()
		{
			const auto div_tests = generate_div_tests_impl();

			size_t largest_remainders_sum = 0;

			// for each div test
			for (const auto& div_test : div_tests)
			{
				// calculate the sum of having every possible bit set, regardless of max_n_of_remainders
				size_t remainders_sum = 0;
				for (size_t i = 0; i < max_pn_bitwidth; ++i)
					remainders_sum += div_test.remainders[i % div_test.n_of_remainders];

				// keep track of the largest sum
				if (remainders_sum > largest_remainders_sum)
					largest_remainders_sum = remainders_sum;
			}

			// + 1 so "lookup[sum]" is always safe.
			return largest_remainders_sum + 1;
		}
		constexpr size_t prime_factor_lookup_size = calculate_prime_factor_lookup_size();

		using prime_lookup_t = narrowest_uint_for_n_bits<div_test::n_of_primes>;
		std::vector<prime_lookup_t> build_prime_factor_lookup()
		{
			std::vector<prime_lookup_t> lookup;
			lookup.reserve(prime_factor_lookup_size);

			// for every possible summation of remainders
			for (size_t i = 0; i < prime_factor_lookup_size; ++i)
			{
				// for each prime
				prime_lookup_t entry = 0;
				for (prime_lookup_t p = 0; p < div_test::n_of_primes; ++p)
				{
					// calculate if that summation of remainders is divisible by that prime
					entry |= (prime_lookup_t((i % small_primes_lookup[p]) == 0) << p);
				}

				lookup.push_back(entry);
			}

			return lookup;
		}

		// Replaces "n % prime[k] == 0" with "lookup[n] & (1 << k)"
		const std::vector<prime_lookup_t> prime_factor_lookup = build_prime_factor_lookup();
	}

	__forceinline bool has_small_prime_factor(const size_t n, const prime_idx_t prime_index)
	{
		return (detail::prime_factor_lookup[n] & (detail::prime_lookup_t(1) << prime_index)) != 0;
	}



	// looping, sorted div tests:

	constexpr size_t div_tests_size = detail::generate_div_tests_impl().size();
	consteval std::array<div_test_t, div_tests_size> generate_div_tests()
	{
		std::array<div_test_t, div_tests_size> div_tests;
		const auto x = detail::generate_div_tests_impl();
		std::copy(x.begin(), x.end(), div_tests.begin());
		return div_tests;
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
