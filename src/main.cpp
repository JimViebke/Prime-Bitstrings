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
	1000000101110111001011011000101000111000010001111 - a p8
	1000000110101000000110000011010001000110011101011 - a p8
	*/
	size_t number = 0b1000000010000011110100010001000101001010110111001;
	static mpz_class mpz_prime; // it's a surprise tool that will help us later
	mpz_prime = 0ull;

	const size_t stopping_point = number + 500'000'000;
	const std::vector<uint8_t> static_sieve = generate_static_sieve();
	/* The number must start on an odd multiple of the sieve size.
	 * To round N to the nearest odd multiple of K:
	 * n -= k;
	 * n -= n % 2k;
	 * n += k; */
	number -= static_sieve.size();
	number -= number % (2 * static_sieve.size());
	number += static_sieve.size();

	// Don't start the clock until here
	auto start = current_time_in_ms();

	for (; number < stopping_point; )
	{
		// perform additional sieving on the static sieve
		std::vector<uint8_t> sieve = static_sieve;
		partial_sieve(number, sieve);

		for (size_t i = 0; i < static_sieve.size(); ++i, number += 2)
		{
			// Bail if this number is already known to have a prime factor < 1000
			if (!sieve[i]) continue;

			// Bail if n does not have a prime number of bits set.
			if ((tiny_primes_lookup() & (1ull << pop_count(number))) == 0) continue;
			// We could use a 64-byte lookup instead of 64-bit lookup:
			// if (!tiny_primes[pop_count(number)]) continue;

			// Bail if gcd(abs(# of even bits - # of odd bits), 1155) is not equal to one.
			const int pca = (int)pop_count(number & 0xAAAAAAAAAAAAAAAA);
			const int pcb = (int)pop_count(number & 0x5555555555555555);
			if (gcd_1155[abs(pca - pcb)] != 1) continue;
			// Instead of "discard if gcd( ... ) != 1", you can "discard the candidate if sa is not a prime greater than 12."
			// These are effectively performance-indentical, as both are some_lookup[abs(pca - pcb)]

			// Bail if n is not prime in base 2
			mpz_prime = number;
			if (!mpir_is_prime(mpz_prime, small_primes_cap)) continue;
			// misof_16k::is_prime(number); // ~2.5x slower
			// misof_262k::is_prime_2_64(number); // ~1.25x slower
			// if (!pk::is_prime(number)) continue; // ~6300x slower (??)
			// if (!pk::fast_is_prime(number)) continue; // still extremely slow

			// convert uint64_t to char array of ['0', '1'...] for MPIR
			char bin_str[64 + 1];
			auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
			*result.ptr = '\0';

			for (int base = 3; ; ++base)
			{
				mpz_prime.set_str(bin_str, base);
				if (!mpir_is_prime(mpz_prime))
				{
					if (base > 8)
					{
						log_result(number, base - 1); // prime up to base - 1, not base
					}

					break;
				}
			}
		}
	}

	std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";
}


int main()
{
	find_multibase_primes();
}
