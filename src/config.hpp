#pragma once

#include <stdint.h>
#include <array>

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
	*/

	const size_t p11 = 0b1000000010000011110100010001000101001010110111001;
	const size_t largest_known_mpb = 0b1000010100100001100110111110000101010011001100111;

	const bool benchmark_mode = true;

	const size_t starting_point = benchmark_mode ? p11 : largest_known_mpb;
	const size_t stopping_point = starting_point + 5'000'000'000;

	// The size of the static sieve is the product of these numbers. Exercise caution.
	constexpr std::array static_sieve_primes{ 3, 5, 7, 11, 13 };

	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(), static_sieve_primes.end(), size_t(1), std::multiplies());

	namespace div_test // trial division tests
	{
		// 2*20 or 2*18 give best performance so far
		const size_t rounds = 2;
		const size_t primes_per_round = 20;
		const size_t n_of_primes = rounds * primes_per_round;

		const size_t up_to_base = 8;
	}

	const size_t smallest_base_to_log = 8;

	const size_t sieve_primes_cap = 1621; // previously 1000
}
