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

	const size_t bm_start = p11; // default: p11
	const size_t bm_size = 50'000'000'000; // default 50'000'000'000
	const size_t bm_stop = bm_start + bm_size;

	constexpr size_t num_threads = 1; // could be passed by command line instead
	constexpr size_t block_size = 1'000'000'000;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull };

	constexpr size_t sieve_alignment = sizeof(__m256i);

	// For 50B benchmark, combinations of 3*43*12907 yield the same pass counts as 1
	constexpr size_t sieve_steps = 3ull * 43ull; // default: 3ull * 43ull

	const size_t sieve_primes_cap = 2200; // default: 2200

	namespace div_test // trial division tests
	{
		constexpr size_t n_of_primes = 32; // default: 32
		constexpr size_t up_to_base = 12; // default: 12

		constexpr size_t max_remainders = 50; // default: 50
	}

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
#else
#define count_passes(...)
#endif

#if 1 // toggle inlining on sieve, popcount, GCD, and div tests
#define tests_are_inlined __forceinline
#else
#define tests_are_inlined __declspec(noinline)
#endif
