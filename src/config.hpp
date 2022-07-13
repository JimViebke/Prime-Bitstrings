#pragma once

#include <array>
#include <numeric>

namespace mbp
{
	/*
	10011110011011110110110011 - p9
	1000000010000011110100010001000101001010110111001 - p11
	*/

	const size_t p11 = 0b1000000010000011110100010001000101001010110111001;

	const bool benchmark_mode = true;

	const size_t bm_start = p11;
	const size_t bm_size = 5'000'000'000;
	const size_t bm_stop = bm_start + bm_size;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull };

	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(), static_sieve_primes.end(), size_t(1), std::multiplies());

	const size_t sieve_primes_cap = 1621; // previously 1000

	namespace div_test // trial division tests
	{
		const size_t n_of_primes = 20;
		const size_t up_to_base = 8;

		// Full div testing should step over the #of primes (starting from 3) with a hardcoded divtest
		constexpr size_t n_of_primes_with_hardcoded_divtests = 3; // 3 == skip 3, 5, 7

		// if up_to_base==8, then (3..n values) == (6 values) == (n + 1 - 3 values)
		constexpr size_t n_of_bases = (up_to_base + 1) - 3;

		// -1 because we skip 2
		constexpr size_t mod_remainders_size = (n_of_primes - n_of_primes_with_hardcoded_divtests - 1) * n_of_bases;

		constexpr size_t max_div_test_remainders = 64; // max 64
	}

	namespace prime_test
	{
		const size_t n_random_bases = 1;
	}

	const char results_path[] = "results.txt";
}
