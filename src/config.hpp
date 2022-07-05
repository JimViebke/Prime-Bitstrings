#pragma once

#include <stdint.h>
#include <array>
#include <numeric>

namespace mbp
{
	/*
	Checkpoints to copy/paste as the starting point:
	10011110011011110110110011 - p9
	1000000010000011110100010001000101001010110111001 - p11
	1000000011110100000000010000000101100110001111111 - a p9
	1000000100010101001111100110110101001001100011101 - a p9
	1000000101000001101101010011100101101110001100111 - a p9
	1000000101000010110111001101110110110000101010011 - a p9
	1000000101000100101111101000110001001111110101001 - a p9
	1000001101100110101110010111111111001111010001011 - a p8
	1000010000010110011000010000011010101110110100001 - a p8
	1000010000100110000000011111000100100001110011101 - a p8
	1000010100100001100110111110000101010011001100111 - a p8
	1000011010100011010110010011101011010011010000011 - a p8
	*/

	const size_t p11 = 0b1000000010000011110100010001000101001010110111001;

	const bool benchmark_mode = false;

	const size_t bm_start = p11;
	const size_t bm_size = 5'000'000'000;
	const size_t bm_stop = bm_start + bm_size;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3ull, 5ull, 7ull, 11ull, 13ull };

	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(), static_sieve_primes.end(), size_t(1), std::multiplies());

	namespace div_test // trial division tests
	{
		const size_t n_of_primes = 20;

		const size_t up_to_base = 8;
	}

	const size_t sieve_primes_cap = 1621; // previously 1000

	const char results_path[] = "results.txt";
}
