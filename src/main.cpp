// Prime Bitstrings.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>
#include <charconv>
#include <bitset>

#include "pk_prime.hpp"

#include "utility.hpp"
#include "math.hpp"
#include "multibase_div_tests.hpp"

void log_time()
{
	time_t timestamp = time(0);
	tm now;
	localtime_s(&now, &timestamp);
	std::cout << std::setfill('0');
	std::cout << std::setw(2) << ((now.tm_hour % 12 == 0) ? 12 : now.tm_hour % 12) << ':' << std::setw(2) << now.tm_min << '\t';
}

void log_result(const mpz_class& n, size_t up_to_base)
{
	log_time();

	std::stringstream ss;
	ss << bin_to_base(n, 10) << " is a p" << up_to_base << " (" << n << ")\n";

	std::cout << ss.str();

	// Only log large primes to file
	if (up_to_base >= 8)
	{
		std::ofstream ofs("results.txt", std::ofstream::app);
		ofs << ss.str();
	}
}

void partial_sieve(const size_t start, std::vector<uint8_t>& sieve)
{
	for (const size_t p : small_primes_lookup)
	{
		// only mark off the small primes not already here
		if (p <= 13) continue;

		// Find out how far past we are from the previous multiple of p.
		// ie 3 == 10 % 7
		// This will always be <= p.
		size_t offset_from_last_pn = (start % p) / 2;
		// Divide by 2 because the sieve only represents odd numbers

		// Now mark false each multiple of p, starting from [0 + p - distance from previous p],
		// where 0 is the n we started with.
		for (size_t i = p - offset_from_last_pn; i < sieve.size(); i += p)
		{
			sieve[i] = false;
		}
	}
}

const std::vector<uint8_t> generate_static_sieve()
{
	const std::vector<size_t> sieve_primes = { 3, 5, 7, 11, 13 };

	std::vector<uint8_t> sieve({ 3 * 5 * 7 * 11 * 13 }, true);

	// for each prime, mark off all multiples
	for (const size_t p : sieve_primes)
		for (size_t i = p; i < sieve.size(); i += p)
			sieve[i] = false;

	return sieve;
}

void calculate_static_sieve_sizes()
{
	size_t sieve_size = 1;
	for (auto p : small_primes_lookup)
	{
		if (p == 2) continue;
		sieve_size *= p;
		std::cout << "Static sieve size would be size " << sieve_size << " using product of primes up to " << p << '\n';
		if (sieve_size > 1'000'000'000) return;
	}
	std::cout << "(no suitable sieve size found)\n";
}

void find_multibase_primes()
{
	gmp_random::r.seed(rand());

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
	*/
	size_t number = 0b1000000010000011110100010001000101001010110111001;
	mpz_class mpz_prime = 0ull; // it's a surprise tool that will help us later

	const size_t stopping_point = number + 5'000'000'000;
	const std::vector<uint8_t> static_sieve = generate_static_sieve();
	std::vector<uint8_t> sieve;
	/* The number must start on an odd multiple of the sieve size.
	 * To round N to the nearest odd multiple of K:
	 * n -= k;
	 * n -= n % 2k;
	 * n += k; */
	number -= static_sieve.size();
	number -= number % (2 * static_sieve.size());
	number += static_sieve.size();

	constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();
	constexpr size_t gcd_1155_lookup = build_gcd_1155_lookup();

	// Dimensions are [base 3..n][primes][residues]
	const std::vector<std::vector<std::vector<uint8_t>>> remainders = generate_remainders_for_bases(12, 40);
	// Dimensions are [base 3..n][bitmasks for p]
	const std::vector<std::vector<size_t>> bitmasks = generate_mod_remainder_bitmasks(remainders);

	// Don't start the clock until here
	const auto start = current_time_in_ms();

	for (; number < stopping_point; )
	{
		// perform additional sieving on the static sieve
		sieve = static_sieve;
		partial_sieve(number, sieve);

		for (size_t i = 0; i < static_sieve.size(); ++i, number += 2)
		{
			// Bail if this number is already known to have a prime factor < 1000
			if (!sieve[i]) continue;

			// Bail if n does not have a prime number of bits set.
			if ((tiny_primes_lookup & (1ull << pop_count(number))) == 0) continue;

			// Bail if gcd(abs(# of even bits - # of odd bits), 1155) is not equal to one.
			const int pca = (int)pop_count(number & 0xAAAAAAAAAAAAAAAA);
			const int pcb = (int)pop_count(number & 0x5555555555555555);
			if ((gcd_1155_lookup & (1ull << abs(pca - pcb))) == 0) continue;

			// Do cheap(er) trial division tests

			const size_t step = 20; // 20*2 or 18*2 give best performance so far
			for (size_t j = 0; j < step * 2; j += step)
			{
				// for each base 3..8
				for (size_t base = 3; base <= 8; ++base)
				{
					// for each small prime
					for (size_t k = j; k < j + step; ++k)
					{
						// mask against bitmask[base][k] to collect residues in each set of positions
						size_t residues = 0;
						for (size_t l = 0; l < remainders[base][k].size(); ++l)
						{
							residues += pop_count(number & (bitmasks[base][k] << l)) * remainders[base][k][l];
						}

						// see if the sum of residues is evenly divisible by a given prime
						if (divides_evenly(residues, k)) { goto done; }
					}
				}
			}

			if (false)
			{
			done:
				continue;
			}

			// Do full primality tests, bail when n is not prime

			if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

			// convert uint64_t to char array of ['0', '1'...] for MPIR
			char bin_str[64 + 1];
			auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
			*result.ptr = '\0';

			mpz_prime.set_str(bin_str, 3);
			if (!mpir_is_prime(mpz_prime)) continue;

			mpz_prime.set_str(bin_str, 4);
			if (!mpir_is_prime(mpz_prime)) continue;

			mpz_prime.set_str(bin_str, 5);
			if (!mpir_is_prime(mpz_prime)) continue;

			mpz_prime.set_str(bin_str, 6);
			if (!mpir_is_prime(mpz_prime)) continue;

			mpz_prime.set_str(bin_str, 7);
			if (!mpir_is_prime(mpz_prime)) continue;

			mpz_prime.set_str(bin_str, 8);
			if (!mpir_is_prime(mpz_prime)) continue;

			mpz_prime.set_str(bin_str, 9);
			if (!mpir_is_prime(mpz_prime)) { log_result(number, 8); continue; }

			mpz_prime.set_str(bin_str, 10);
			if (!mpir_is_prime(mpz_prime)) { log_result(number, 9); continue; }

			mpz_prime.set_str(bin_str, 11);
			if (!mpir_is_prime(mpz_prime)) { log_result(number, 10); continue; }

			mpz_prime.set_str(bin_str, 12);
			if (!mpir_is_prime(mpz_prime)) { log_result(number, 11); continue; }

			mpz_prime.set_str(bin_str, 13);
			if (!mpir_is_prime(mpz_prime)) { log_result(number, 12); continue; }
		}
	}

	std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";
}


int main()
{
	find_multibase_primes();
}
