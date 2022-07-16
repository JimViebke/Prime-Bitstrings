#pragma once

#include <vector>

#include "pk_prime.hpp"
#include "math.hpp"

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
#if display_unused_div_tests
		bool used = false;
		base_t base;
		uint32_t hits = 0;
#endif
		prime_idx_t prime_idx = 0;
		n_of_remainders_t n_of_remainders = 0; // is also the index of the req'd bitmask
		std::array<remainder_t, max_remainders + 1> remainders{ 0 };
	};

	namespace detail
	{
		constexpr std::vector<div_test_t> generate_div_tests_impl()
		{
			std::vector<div_test_t> div_tests;

			for (size_t i = n_of_primes_with_hardcoded_divtests + 1; i < n_of_primes; ++i) // for each small prime
			{
				for (size_t base = 3; base <= up_to_base; ++base) // for each base 3..n
				{
#if display_unused_div_tests
					div_test_t dt{ .base = base_t(base), .prime_idx = prime_idx_t(i) };
#else
					const auto p = small_primes_lookup[i];

					if (base == 4 && p == 3) continue; // base 4^n % 3 unused
					if (base == 5 && p == 3) continue; // base 5^n % 3 unused
					if (base == 7 && p == 3) continue; // base 7^n % 3 unused
					if (base == 8 && p == 3) continue; // base 8^n % 3 unused
					if (base == 10 && p == 3) continue; // base 10^n % 3 unused
					if (base == 11 && p == 3) continue; // base 11^n % 3 unused

					if (base == 4 && p == 5) continue; // base 4^n % 5 unused
					if (base == 6 && p == 5) continue; // base 6^n % 5 unused
					if (base == 7 && p == 5) continue; // base 7^n % 5 unused
					if (base == 9 && p == 5) continue; // base 9^n % 5 unused
					if (base == 11 && p == 5) continue; // base 11^n % 5 unused

					if (base == 6 && p == 7) continue; // base 6^n % 7 unused
					if (base == 8 && p == 7) continue; // base 8^n % 7 unused
					if (base == 9 && p == 7) continue; // base 9^n % 7 unused

					if (base == 10 && p == 11) continue; // base 10^n % 11 unused

					// If two div tests are identical, remove one
					if (base == 8 && p == 5) continue; //  base 3^n % 5 is congruent to  8^n % 5
					if (base == 10 && p == 7) continue; // base 3^n % 7 is congruent to 10^n % 7
					if (base == 11 && p == 7) continue; // base 4^n % 7 is congruent to 11^n % 7


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
						  //return a.n_of_remainders < b.n_of_remainders;
						  //return a.prime_idx < b.prime_idx;

						  return
							  size_t(a.n_of_remainders) * small_primes_lookup[a.prime_idx] <
							  size_t(b.n_of_remainders) * small_primes_lookup[b.prime_idx];
					  });

			//for (auto& dt : div_tests)
			//	std::cout << "b" << size_t(dt.base) << " % " << small_primes_lookup[dt.prime_idx] << " \t" << size_t(dt.n_of_terms) << " terms\n";

			return div_tests;
		}

		// bitmask for all base^n mod prime
		template<size_t base, size_t prime>
		struct bitmask_for
		{
			static constexpr size_t val = []
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
			}();
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
			static constexpr size_t val = []
			{
				for (size_t i = 1; i < 64; ++i)
					if ((bitmask >> i) & 1)
						return i;
			}();
		};

		template<size_t prime>
		struct get_prime_index
		{
			static constexpr size_t idx = []
			{
				for (size_t i = 0; i < 64; ++i)
					if (small_primes_lookup[i] == prime)
						return i;
			}();
		};

		constexpr size_t calculate_prime_factor_lookup_size()
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

	// Suppress warnings about bitmasks having upper bits moved
#pragma warning(push)
#pragma warning(disable: 26450)

	// unrolled div tests, templated version:

	template<size_t base>
	struct in_base
	{
		static constexpr size_t val = base;
		static constexpr bool is(size_t b) { return b == base; }
	};

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 5 && base == 3)
		/*__forceinline*/ bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 4);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;

		return has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 7 && base == 3)
		/*__forceinline*/ bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<base, 4, divisor>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<base, 5, divisor>::rem;

		return has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 7 && base == 4)
		/*__forceinline*/ bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 3);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;

		return has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 7 && base == 5)
		/*__forceinline*/ bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<base, 4, divisor>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<base, 5, divisor>::rem;

		return has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	//template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 11 && base == 8)
	//	__forceinline bool is_divisible_by(const size_t number)
	//{
	//	using namespace detail;
	//	constexpr size_t mask = bitmask_for<base, divisor>::val;
	//	static_assert(period_of<mask>::val == 10);
	//	size_t rem = 0;
	//	rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
	//	rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
	//	rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
	//	rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;
	//	rem += pop_count(number & (mask << 4)) * pow_mod<base, 4, divisor>::rem;
	//	rem += pop_count(number & (mask << 5)) * pow_mod<base, 5, divisor>::rem;
	//	rem += pop_count(number & (mask << 6)) * pow_mod<base, 6, divisor>::rem;
	//	rem += pop_count(number & (mask << 7)) * pow_mod<base, 7, divisor>::rem;
	//	rem += pop_count(number & (mask << 8)) * pow_mod<base, 8, divisor>::rem;
	//	rem += pop_count(number & (mask << 9)) * pow_mod<base, 9, divisor>::rem;

	//	return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	//}

	// unrolled div tests, recursive templated version:

	template<size_t divisor, size_t base, size_t mask, size_t place_value>
	__forceinline void recursive_is_divisible_by(size_t& rem, const size_t number)
	{
		rem += pop_count(number & (mask << place_value)) * detail::pow_mod<base, place_value, divisor>::rem;
		if constexpr (place_value > 0)
			recursive_is_divisible_by<divisor, base, mask, place_value - 1>(rem, number);
		return;
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val>
	__forceinline bool recursive_is_divisible_by(const size_t number)
	{
		static_assert(base >= 3 && base <= up_to_base);
		constexpr size_t bitmask = detail::bitmask_for<base, divisor>::val;
		size_t rem = 0;
		recursive_is_divisible_by<divisor, base, bitmask, detail::period_of<bitmask>::val - 1>(rem, number);
		return has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

#pragma warning(pop)

	// looping div tests, but sorted:

	constexpr size_t div_tests_size = detail::generate_div_tests_impl().size();
	constexpr std::array<div_test_t, div_tests_size> generate_div_tests()
	{
		std::array<div_test_t, div_tests_size> div_tests;
		const auto x = detail::generate_div_tests_impl();
		std::copy(x.begin(), x.end(), div_tests.begin());
		return div_tests;
	}
}
