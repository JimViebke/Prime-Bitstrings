#pragma once

#include <array>
#include <numeric>

#define display_unused_div_tests 0

namespace mbp
{
	/*
	10011110011011110110110011 - p9
	1000000010000011110100010001000101001010110111001 - p11
	*/

	const size_t p11 = 0b1000000010000011110100010001000101001010110111001;

	const bool benchmark_mode = false;

	const size_t bm_start = p11;
	const size_t bm_size = 5'000'000'000;
	const size_t bm_stop = bm_start + bm_size;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull };

	const size_t sieve_primes_cap = 2200; // 2200

	namespace div_test // trial division tests
	{
		const size_t n_of_primes = 32; // 32
		const size_t up_to_base = 11; // 11

		constexpr size_t max_remainders = 42; // 42

		// Full div testing should step over the #of primes (starting from 3) with a hardcoded divtest
		constexpr size_t n_of_primes_with_hardcoded_divtests = 0; // 3 == skip 3, 5, 7

		// if up_to_base==8, then (3..n values) == (6 values) == (n + 1 - 3 values)
		constexpr size_t n_of_bases = (up_to_base + 1) - 3;

		// Probably don't touch
		constexpr size_t max_pn_bitwidth = 50;
		static_assert(max_pn_bitwidth > std::bit_width(p11) &&
					  max_pn_bitwidth <= std::numeric_limits<size_t>::digits);
	}

	namespace prime_test
	{
		const size_t n_random_bases = 1;
	}

	const char results_path[] = "results.txt";
}
