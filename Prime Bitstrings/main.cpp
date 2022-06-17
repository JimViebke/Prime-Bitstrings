// Prime Bitstrings.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>
#include <charconv>

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

	std::ofstream ofs("../results.txt", std::ofstream::app);
	ofs << ss.str();
}

void partial_sieve(const size_t& start, std::vector<uint8_t>& sieve)
{
	for (const size_t p : small_primes_lookup)
	{
		// Find out how far past we are from the previous multiple of p.
		// ie 3 == 10 % 7
		// This will always be <= p.
		size_t offset_from_last_pn = start % p;

		// Now mark false each multiple of p, starting from [0 + p - distance from previous p],
		// where 0 is the n we started with.
		for (size_t i = p - offset_from_last_pn; i < sieve.size(); i += p)
		{
			sieve[i] = false;
		}
	}
}

void find_multibase_primes()
{
	gmp_random::r.seed(rand());

	auto start = current_time_in_ms();

	/*
	Checkpoints to copy/paste as the starting point:
	10011110011011110110110011 - p9
	1000000010000011110100010001000101001010110111001 - p11
	1000000011110100000000010000000101100110001111111 - a p9
	1000000100010101001111100110110101001001100011101 - a p9
	1000000101000001101101010011100101101110001100111 - a p9
	1000000101000010110111001101110110110000101010011 - a p9
	1000000101000100101111101000110001001111110101001 - p9
	1000000101110111001011011000101000111000010001111 - p8
	1000000101110111001011011000101000111000010001111
	*/
	size_t number = 0b1000000010000011110100010001000101001010110111001;

	const size_t stopping_point = number + 500'000'000;
	const size_t sieve_size = 1'000'000;

	for (; number < stopping_point; )
	{
		std::vector<uint8_t> sieve(sieve_size, true);
		partial_sieve(number, sieve);

		for (size_t i = 0; i < sieve_size; ++i, ++number)
		{
			// Bail if this number is already known to have a factor < 1000
			if (!sieve[i]) continue;

			// Bail if n does not have a prime number of bits set.
			if ((tiny_primes_lookup() & (1ull << pop_count(number))) == 0) continue;
			// if we were to have tiny_primes be a 64-byte lookup instead of a 64-bit number:
			// if (!tiny_primes[pop_count(number)]) continue;

			// Bail if gcd(abs(# of even bits - # of odd bits), 1155) is not equal to one.
			// Interestingly, this prevents us from detecting p7s, but not p8s+.
			const int pca = (int)pop_count(number & 0xAAAAAAAAAAAAAAAA);
			const int pcb = (int)pop_count(number & 0x5555555555555555);
			if (gcd_1155[abs(pca - pcb)] != 1) continue;

			// Instead of "discard if gcd( ... ) != 1", you can "discard the candidate if sa is not a prime greater than 12." 
			// These are effectively performance-indentical, as both involve some_lookup[ abs(pca - pcb) ].

			// Bail if n is not prime in this base (base 2)
			if (!mpir_is_prime(number, small_primes_cap)) continue; // can we do this natively without calling the lib?
			// mpir_is_prime(number, small_primes_cap); // ~20 s
			// misof_16k::is_prime(number); // ~50 s
			// misof_262k::is_prime_2_64(number); // ~24 s

			// convert uint64_t to char array of ['0', '1'...] for mpz_class
			char bin_str[64 + 1];
			auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
			*result.ptr = '\0';

			for (int base = 3; ; ++base)
			{
				const auto number_to_base = mpz_class{ bin_str, base };
				if (!mpir_is_prime(number_to_base, small_primes_cap))
				{
					// We just failed on "base", so "base > n" only logs
					// results that reached equal to n or higher
					if (base > 8)
					{
						log_result(number, base - 1);
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

	// find_p2_8();
}
