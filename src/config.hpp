#pragma once

#include <array>
#include <numeric>

#include "util/utility.hpp"

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
	p12		10100001011000101000110101011011011101111110100101011
	*/

	const size_t p11 = 0b1000000010000011110100010001000101001010110111001;
	const size_t p12 = 0b10100001011000101000110101011011011101111110100101011;

	const bool benchmark_mode = false;

	const size_t bm_start = p11; // default: p11
	const size_t bm_size = 50'000'000'000; // default 50'000'000'000
	const size_t bm_stop = bm_start + bm_size;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull, 17ull };

	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(),
														 static_sieve_primes.end(), 1ull, std::multiplies());

	constexpr size_t sieve_steps = 2ull; // default 2

	const size_t sieve_primes_cap = 937; // default: 937

	// trial division tests in bases 3..n
	namespace div_test
	{
		constexpr size_t n_of_primes = 32; // default: 32
		constexpr size_t up_to_base = 13; // default: 13

		constexpr size_t reorder_interval = 10'000'000'000; // default: 10 B
	}

	// full primality testing in bases 3..n
	namespace prime_test
	{
		const size_t n_random_bases = 1;
	}

	const char results_filename[] = "results.txt";
}

#define analyze_div_tests 0
#define suppress_extra_div_tests 1



#ifdef __INTELLISENSE__
#define use_constexpr
#define use_consteval
#else
#define use_constexpr constexpr
#define use_consteval consteval
#endif

#if analyze_div_tests
#define div_test_const
#define div_test_constexpr
#else
#define div_test_const const
#define div_test_constexpr use_constexpr
#endif

#if 0 // accumulate and print pass counts
#define count_passes(...) __VA_ARGS__
#define log_pass_counts(str, count, count_before) { \
	auto w = std::setw; \
	std::cout << str << w(10) << count << " (removed ~" << w(3) << 100 - (count * 100 / count_before) << "%) (~" << w(6) << count / passes << " candidates pass per main loop iteration)\n"; \
	pc_hash = util::hash(pc_hash ^ count); }
#else
#define count_passes(...)
#define log_pass_counts(...)
#endif

#if 1 // toggle inlining on sieve, GCD, and div tests
#define tests_are_inlined __forceinline
#else
#define tests_are_inlined __declspec(noinline)
#endif
