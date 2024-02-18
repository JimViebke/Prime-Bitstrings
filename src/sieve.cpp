
#include <iostream>

#include "math/math.hpp"
#include "sieve.hpp"

namespace mbp::prime_sieve
{
	std::array<sieve_offset_t, n_of_vector_sieve_primes> sieve_offsets_cache;

	void set_up_sieve_offsets_cache(const uint64_t start)
	{
		static_assert(sizeof(sieve_offset_t) >= sizeof(sieve_prime_t));

		// start with the first prime not in the static sieve
		auto prime_it = small_primes_lookup.begin() + static_sieve_primes.size() + 1;
		for (size_t i = 0; i < sieve_offsets_cache.size(); ++i)
		{
			const size_t prime = *prime_it++;
			size_t stride = prime;

			if (prime > largest_aligned_vector_sieve_prime &&
				prime <= largest_vector_sieve_prime)
			{
				stride *= 15;
			}

			// find the distance to the next multiple of the stride
			size_t n = stride - (start % stride);

			// make sure start + n is an odd multiple of the stride (start is always odd)
			if (n % 2 == 1)
				n += stride;

			if (n == 2 * stride) // handle edge cases where start % stride == 0
				n = 0;

			n /= 2; // convert from distance to index

			if (prime > largest_aligned_vector_sieve_prime &&
				prime <= largest_vector_sieve_prime)
			{
				// We sieve by strides of 8*15*p, starting with a bit offset of 0.
				// Advance by 15*p until we have this alignment.
				while (n % 8 != 0)
					n += stride;

				// if we've reached the second multiple of 8*15*p, step back to the first
				n = util::min(n, n - 8ull * stride);
			}

			sieve_offsets_cache[i] = sieve_offset_t(n);
		}
	}

	void dump_sieve_offsets_cache()
	{
		std::cout << "Dumping sieve_offsets_cache:\n";

		for (size_t i = 0; i < sieve_offsets_cache.size(); ++i)
		{
			std::cout << i << '\t' << sieve_offsets_cache[i] << '\n';
		}
	}

	void detail::verify_sieve_offset_cache(const uint64_t start)
	{
		for (size_t i = static_sieve_primes.size() + 1; small_primes_lookup[i] <= largest_sieve_prime; ++i)
		{
			const uint64_t prime = small_primes_lookup[i];
			const size_t offset = sieve_offsets_cache[i];

			if (prime >= 17 && prime <= largest_aligned_vector_sieve_prime) // handle separately
			{
				if (offset >= prime)
				{
					std::cout << "offset == " << offset << ", should be less than prime (" << prime << "), start = " << start << '\n';
					dump_sieve_offsets_cache();
					std::cin.ignore();
				}

				if ((start + (2 * offset)) % prime != 0)
				{
					std::cout << "(start + (2 * offset)) % prime == " << (start + (2 * offset)) % prime << ", should be 0 (offset == " << offset << "), prime = " << prime << ", start = " << start << '\n';
					dump_sieve_offsets_cache();
					std::cin.ignore();
				}

				continue;
			}

			const uint64_t p15 = 15ull * prime;

			// 1. start + (2 * offset) should be evenly divisible by 15*p
			if ((start + (2 * offset)) % p15 != 0)
			{
				std::cout << "(start + (2 * offset)) % p15 == " << (start + (2 * offset)) % p15 << ", should be 0\n";
				std::cout << "start == " << start << '\n';
				std::cout << "offset == " << offset << '\n';
				dump_sieve_offsets_cache();
				std::cin.ignore();
			}

			// 2. For vectorized sieving, the offset should be evenly divisible by 8.
			if (prime <= largest_vector_sieve_prime && offset % 8 != 0)
			{
				std::cout << "offset % 8 == " << offset % 8 << ", should be 0\n";
				dump_sieve_offsets_cache();
				std::cin.ignore();
			}

			// 3. The offset should be less than 8 * 15*p
			if (offset >= 8 * p15)
			{
				std::cout << "prime == " << prime << '\n';
				std::cout << "15*p == " << p15 << '\n';
				std::cout << "8 * 15*p == " << 8 * p15 << '\n';
				std::cout << "offset == " << offset << ", should be < " << 8 * p15 << '\n';
				dump_sieve_offsets_cache();
				std::cin.ignore();
			}
		}
	}

}
