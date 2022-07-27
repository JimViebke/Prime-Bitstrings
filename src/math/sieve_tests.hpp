#pragma once

#include <iostream>
#include <vector>

#include "franken_mpir.hpp"

namespace mbp
{
	void sieve_of_euler_2()
	{
		// An optimization of the basic sieve of Eratosthenes.
		// Example:
		// When eliminating products of 7, don't mark off [14, 21, 28, ...]
		// Instead, mark off [7*7, 7*11, 7*13, ...]
		// This eliminates each composite number only once

		std::vector<uint8_t> sieve(200, true);
		sieve[0] = sieve[1] = false;

		std::vector<size_t> primes;
		primes.reserve(sieve.size());

		// for i through N
		for (size_t i = 2; i < sieve.size(); ++i)
		{
			if (sieve[i])
				primes.push_back(i);

			// for each known prime, until we find a prime that divides i:
			for (size_t j = 0; j < primes.size(); ++j)
			{
				if (primes[j] * i >= sieve.size()) break;

				sieve[primes[j] * i] = false;
				std::cout << "Marking off " << i << '*' << primes[j] << " (" << i * primes[j] << ")\n";

				if (i % primes[j] == 0) break;
			}

			std::cout << '\n';
		}

		// std::cout << "The following integers have no prime factors below " << small_primes_lookup.back() << ':';

		for (const auto p : primes)
		{
			std::cout << p << ' ';
		}
	}

	void sieve_tests_3()
	{
		std::vector<uint8_t> sieve(100'000'000, true);
		sieve[0] = sieve[1] = false;
		const size_t stop = size_t(franken_boost::sqrt(sieve.size())) + 1;

		size_t steps = 0;
		size_t primes = 0;

		for (const size_t p : small_primes_lookup)
		{
			if (p > stop) break;

			for (size_t i = p + p; i < sieve.size(); i += p)
			{
				sieve[i] = false;
				++steps;
			}
		}

		for (size_t i = 0; i < sieve.size(); ++i)
			if (sieve[i])
				++primes;

		std::cout << primes << " primes, " << steps << " steps\n";
	}

	void sieve_tests_3_1()
	{
		std::vector<uint8_t> sieve(100'000'000, true);
		sieve[0] = sieve[1] = false;
		const size_t stop = size_t(franken_boost::sqrt(sieve.size())) + 1;

		size_t steps = 0;
		size_t primes = 0;

		for (const size_t p : small_primes_lookup)
		{
			if (p > stop) break;

			for (size_t i = p * p; i < sieve.size(); i += p)
			{
				sieve[i] = false;
				++steps;
			}
		}

		for (size_t i = 0; i < sieve.size(); ++i)
			if (sieve[i])
				++primes;

		std::cout << primes << " primes, " << steps << " steps\n";
	}
}
