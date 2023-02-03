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

	template<bool on_fast_path>
	size_t* branchless_div_tests(size_t* const candidates_begin,
								 size_t* const candidates_end,
								 const size_t n_of_tests);

	template<bool on_fast_path>
	size_t* branching_div_tests(size_t* input,
								const size_t* const candidates_end,
								const size_t start_offset);

	void update_div_test_order();
	void print_div_tests();
	void run_div_test_analysis(const size_t number);

	namespace detail
	{
		// consteval std::vector<div_test_t> generate_div_tests_impl();

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

	namespace detail
	{
		using prime_lookup_t = util::narrowest_uint_for_n_bits<div_test::n_of_primes>;

		std::vector<prime_lookup_t> build_prime_factor_lookup_old();
		std::vector<prime_lookup_t> build_prime_factor_lookup_new();

		const std::vector<prime_lookup_t> prime_factor_lookup = build_prime_factor_lookup_old();

		extern std::array<std::vector<uint8_t>, div_test::n_of_primes> indivisible_by;
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
		constexpr __forceinline void get_sum_of_rems(size_t& rem, const size_t number)
		{
			rem += pop_count(number & (mask << place_value)) * pow_mod<base, place_value, divisor>::rem;
			if constexpr (place_value > 0)
				get_sum_of_rems<divisor, base, mask, place_value - 1>(rem, number);
		}

		template<size_t divisor, typename base_t>
		constexpr __forceinline size_t get_sum_of_rems(const size_t number)
		{
			constexpr size_t base = base_t::val;
			static_assert(base >= 3 && base <= up_to_base);
			constexpr size_t bitmask = bitmask_for<base, divisor>::val;
			size_t rem = 0;
			get_sum_of_rems<divisor, base, bitmask, period_of<bitmask>::val - 1>(rem, number);
			return rem;
		}
	}

	template<size_t divisor, typename base_t>
	__forceinline size_t get_upper_sum_of_rems(const size_t number)
	{
		constexpr size_t upper_bits_mask = size_t(-1) << 32;
		return detail::get_sum_of_rems<divisor, base_t>(number & upper_bits_mask);
	}

	template<size_t divisor, typename base_t>
	__forceinline bool recursive_is_divisible_by(const size_t number)
	{
		size_t rem = detail::get_sum_of_rems<divisor, base_t>(number);
		return has_small_prime_factor(rem, get_prime_index<divisor>::idx);
	}

#pragma warning(pop)
}
