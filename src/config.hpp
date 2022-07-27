#pragma once

#include <array>
#include <numeric>

namespace mbp
{
	/*
	p2		10
	p3		10
	p4		101111
	p5		10010111
	p6		111110100001
	p7		111110100001
	p8		11000011101101111
	p9		10011110011011110110110011
	p10		110100000010101111110001010011001110001
	p11		1000000010000011110100010001000101001010110111001
	p12
	*/

	const size_t p11 = 0b1000000010000011110100010001000101001010110111001;

	const bool benchmark_mode = false;

	const size_t bm_size = 5'000'000'000;
	const size_t bm_start = p11; // default: p11
	const size_t bm_stop = bm_start + bm_size;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull };

	// The size of the factorization wheel sieve is the product of these numbers. Exercise caution.
	constexpr std::array wheel_primes{ 2ull, 3ull, 5ull, 7ull, 11ull, 13ull };

	const size_t sieve_primes_cap = 2200; // default: 2200

	namespace div_test // trial division tests
	{
		const size_t n_of_primes = 32; // default: 32
		const size_t up_to_base = 11; // default: 11

		constexpr size_t max_remainders = 42; // default: 42

		// Probably don't touch
		constexpr size_t max_pn_bitwidth = 50;
		static_assert(max_pn_bitwidth > std::bit_width(p11) &&
					  max_pn_bitwidth <= std::numeric_limits<size_t>::digits);
	}

	namespace prime_test
	{
		const size_t n_random_bases = 1;
	}

	const char results_filename[] = "results.txt";
}

#define analyze_div_tests 0
#define suppress_extra_div_tests 1

#define USE_UNCACHED 1 // cache popcounts?

#if analyze_div_tests
#define use_constexpr
#else
#define use_constexpr constexpr
#endif
