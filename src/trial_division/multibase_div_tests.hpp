#pragma once

/*
Determine if a bitstring has a small prime divisor in base b.

Given any bitstring, the value of each digit is either base^(digit position) or 0. For example, in base 7, bits 0, 1, and 2 can only represent 7^0, 7^1, 7^2... or 0.

Therefore, for base b, place value n, and a small prime p, we calculate and store b^n % p. These remainders eventually follow a repeating pattern, often shorter than the chosen n, so we only need to store k unique remainders. At runtime, for each of the k remainders, use a bitmask to select every k-th bit, perform a popcount, and multiply the count by that k-th remainder. Adding these k results together gives us a small sum of remainders, s. If s is evenly divisible by p, then the bitstring is divisible by p in base b. Discard it.
*/

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "../config.hpp"
#include "../math/math.hpp"
#include "../math/pk_prime.hpp"
#include "../util/utility.hpp"
#include "types.hpp"

namespace mbp::div_test
{
	using div_tests_t = std::vector<div_test::div_test_t>;

	class full_div_tests
	{
	public:
		full_div_tests();

		template<bool on_fast_path>
		uint64_t* branchless_div_tests(
			uint64_t* const candidates_begin,
			uint64_t* const candidates_end,
			const size_t n_of_tests);

		template<bool on_fast_path>
		uint64_t* branching_div_tests(
			uint64_t* input,
			const uint64_t* const candidates_end,
			const size_t start_offset);

		void update_div_test_order();
		void print_div_tests();
		void run_div_test_analysis(const uint64_t number);

	private:
		void permute_div_tests();

		div_tests_t div_tests;

		std::vector<std::array<remainder_t, 64>> permuted_div_tests;
	};



	namespace detail
	{
		// Bitmask for any base^n % prime.
		template<size_t base, size_t prime>
		constexpr uint64_t bitmask_for{ []() consteval {
			// Calculate the period of base^n % prime.
			size_t i = 0;
			for (; i < 64; ++i)
			{
				remainder_t rem = remainder_t(pk::powMod(base, i, prime));
				// Break when the pattern repeats.
				if (rem == 1 && i > 0) break;
			}

			// Mark every i-th bit starting from bit 0.
			uint64_t bitmask = 0;
			for (size_t j = 0; j < 64 && i < 64; j += i)
			{
				bitmask <<= i;
				bitmask |= 1;
			}

			return bitmask;
		}() };

		template<uint64_t bitmask>
		constexpr size_t period_of{ []() consteval {
			for (size_t i = 1; i < 64; ++i)
				if ((bitmask >> i) & 1)
					return i;
			std::unreachable();
		}() };

		template<size_t prime>
		constexpr size_t prime_index_of{ []() consteval {
			for (size_t i = 0; i < 64; ++i)
				if (small_primes_lookup[i] == prime)
					return i;
			std::unreachable();
		}() };
	}

	namespace detail
	{
		using prime_lookup_t = util::narrowest_uint_for_n_bits<div_test::n_of_primes>;

		std::vector<prime_lookup_t> build_prime_factor_lookup_old();
		std::vector<prime_lookup_t> build_prime_factor_lookup_new();

		const std::vector<prime_lookup_t> prime_factor_lookup = build_prime_factor_lookup_old();

		extern const std::array<std::vector<uint8_t>, div_test::n_of_primes> indivisible_by;
	}

	// Replaces "n % prime[idx] == 0" with "lookup[n] & (1 << idx)", usually as bittest + cmov
	__forceinline bool has_small_prime_factor(const size_t n, const prime_idx_t prime_index)
	{
		return (detail::prime_factor_lookup[n] &
			(detail::prime_lookup_t(1) << prime_index)) != 0;
	}
}
