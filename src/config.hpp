#pragma once

#include <array>
#include <functional>
#include <numeric>
#include <string_view>

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

	inline constexpr size_t p11{ 0b1000000010000011110100010001000101001010110111001 };
	inline constexpr size_t p12{ 0b10100001011000101000110101011011011101111110100101011 };

	inline constexpr size_t bm_start{ p11 };
	inline constexpr size_t bm_size{ 500'000'000'000 };
	inline constexpr size_t bm_stop{ bm_start + bm_size };

	namespace prime_sieve
	{
		inline constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull };
		inline constexpr size_t product_of_static_sieve_primes{ std::accumulate(
				static_sieve_primes.begin(),
				static_sieve_primes.end(),
				1ull, std::multiplies()) };
		inline constexpr size_t static_sieve_size{ 2 * 8 * product_of_static_sieve_primes };

		inline constexpr size_t steps{ 16 };
		inline constexpr size_t largest_aligned_vector_sieve_prime{ 47 };
		inline constexpr size_t largest_vector_sieve_prime{ 79 };
		inline constexpr size_t largest_sieve_prime{ 263 };

		inline constexpr double vector_density_threshold{ 0.0035 };
		inline constexpr double unaligned_vector_density_threshold{ 0.006 };
		inline constexpr double scalar_density_threshold{ 0.030 };
	}

	// trial division tests in bases 2..n
	namespace div_test
	{
		inline constexpr size_t n_of_primes{ 54 };
		inline constexpr size_t up_to_base{ 13 };

		inline constexpr size_t n_of_branchless_tests{ 75 };

		inline constexpr size_t reorder_interval{ 10'000'000'000 };
	}

	// full primality testing in bases 3..n
	namespace prime_test
	{
		inline constexpr size_t mpir_trial_div_cap{ 937 };
		inline constexpr size_t n_random_bases{ 1 };
	}

	inline constexpr std::string_view results_filename{ "results.txt" };
}

#define analyze_div_tests 0
#define suppress_extra_div_tests 1



#if 0
// Accumulate and print pass counts.
#include <iomanip>
#include <iostream>
#define count_passes(...) __VA_ARGS__
#define log_pass_counts(str, count, count_before) { \
	auto w = std::setw; \
	std::cout << str << w(10) << count << " (removed ~" << w(3) << 100 - (count * 100 / count_before) << "%) (~" << w(6) << count / passes << " candidates pass per main loop iteration)\n"; \
	pc_hash = util::hash(pc_hash ^ count); }
#else
#define count_passes(...)
#define log_pass_counts(...)
#endif

#if 1
#define inline_toggle __forceinline
#else
#define inline_toggle __declspec(noinline)
#endif
