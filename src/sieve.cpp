
#include <iostream>

#include "math/math.hpp"
#include "sieve.hpp"

namespace mbp::prime_sieve
{
	const sieve_container generate_static_sieve()
	{
		sieve_container sieve{};
		sieve.set_all();

		// for each prime, mark off all multiples
		for (const auto p : static_sieve_primes)
			for (size_t i = 0; i < sieve.size(); i += p)
				sieve.clear_bit(i);

		return sieve;
	}

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



	template<size_t p>
	__forceinline void clear_two_adjacent(sieve_container& sieve, size_t offset)
	{
		// Mark two (adjacent) multiples of p, using one scalar write if possible

		if constexpr (p + 7 < 64)
		{
			constexpr uint64_t mask = ((1ull << p) | 1ull);

			uint8_t* out = sieve.data() + offset / 8;

			if constexpr (p + 7 < 32)
			{
				*(uint32_t*)out &= ~(mask << (offset % 8));
			}
			else
			{
				*(uint64_t*)out &= ~(mask << (offset % 8));
			}
		}
		else
		{
			// two writes are required; use 8-bit writes
			sieve.clear_bit(offset);
			sieve.clear_bit(offset + p);
		}
	}

	// given a ptr to the sieve byte containing n*15p + 1p, where (n*15p + 1p) % 8 == 0:
	//    mark off (n+m)*15p + m1*p for m in range 0..7
	template<size_t p, size_t m1>
		requires (m1 < 15)
	__forceinline void clear_next_eight(uint8_t* ptr_to_1p)
	{
		// how many bytes from sieve_ptr is j + m1*p?
		constexpr size_t byte_offset_0 = ((m1 - 1) * p + (0 * 15 * p)) / 8;
		constexpr size_t byte_offset_1 = ((m1 - 1) * p + (1 * 15 * p)) / 8;
		constexpr size_t byte_offset_2 = ((m1 - 1) * p + (2 * 15 * p)) / 8;
		constexpr size_t byte_offset_3 = ((m1 - 1) * p + (3 * 15 * p)) / 8;
		constexpr size_t byte_offset_4 = ((m1 - 1) * p + (4 * 15 * p)) / 8;
		constexpr size_t byte_offset_5 = ((m1 - 1) * p + (5 * 15 * p)) / 8;
		constexpr size_t byte_offset_6 = ((m1 - 1) * p + (6 * 15 * p)) / 8;
		constexpr size_t byte_offset_7 = ((m1 - 1) * p + (7 * 15 * p)) / 8;

		// how many bits into the above byte is j + m1*p?
		constexpr size_t bit_offset_0 = ((m1 - 1) * p + (0 * 15 * p)) % 8;
		constexpr size_t bit_offset_1 = ((m1 - 1) * p + (1 * 15 * p)) % 8;
		constexpr size_t bit_offset_2 = ((m1 - 1) * p + (2 * 15 * p)) % 8;
		constexpr size_t bit_offset_3 = ((m1 - 1) * p + (3 * 15 * p)) % 8;
		constexpr size_t bit_offset_4 = ((m1 - 1) * p + (4 * 15 * p)) % 8;
		constexpr size_t bit_offset_5 = ((m1 - 1) * p + (5 * 15 * p)) % 8;
		constexpr size_t bit_offset_6 = ((m1 - 1) * p + (6 * 15 * p)) % 8;
		constexpr size_t bit_offset_7 = ((m1 - 1) * p + (7 * 15 * p)) % 8;

		ptr_to_1p[byte_offset_0] &= ~(1ull << bit_offset_0);
		ptr_to_1p[byte_offset_1] &= ~(1ull << bit_offset_1);
		ptr_to_1p[byte_offset_2] &= ~(1ull << bit_offset_2);
		ptr_to_1p[byte_offset_3] &= ~(1ull << bit_offset_3);
		ptr_to_1p[byte_offset_4] &= ~(1ull << bit_offset_4);
		ptr_to_1p[byte_offset_5] &= ~(1ull << bit_offset_5);
		ptr_to_1p[byte_offset_6] &= ~(1ull << bit_offset_6);
		ptr_to_1p[byte_offset_7] &= ~(1ull << bit_offset_7);
	}

	// given a ptr to the sieve byte containing n*15p + 1p, where (n*15p + 1p) % 8 == 0:
	//    mark off (n+m)*15p + m1*p and (n+m)*15p + m2*p for m in range 0..7
	template<size_t p, size_t m1, size_t m2>
		requires (m1 + 1 == m2) && (m2 < 15)
	__forceinline void clear_next_eight(uint8_t* ptr_to_1p)
	{
		// Mark eight pairs of adjacent multiples of p, using one scalar write if possible

		if constexpr (p + 7 < 64)
		{
			using write_t = util::narrowest_uint_for_n_bits<p + 7>;

			constexpr uint64_t mask = ((1ull << p) | 1ull);

			// how many bytes from sieve_ptr is j + m1*p?
			constexpr size_t byte_offset_0 = ((m1 - 1) * p + (0 * 15 * p)) / 8;
			constexpr size_t byte_offset_1 = ((m1 - 1) * p + (1 * 15 * p)) / 8;
			constexpr size_t byte_offset_2 = ((m1 - 1) * p + (2 * 15 * p)) / 8;
			constexpr size_t byte_offset_3 = ((m1 - 1) * p + (3 * 15 * p)) / 8;
			constexpr size_t byte_offset_4 = ((m1 - 1) * p + (4 * 15 * p)) / 8;
			constexpr size_t byte_offset_5 = ((m1 - 1) * p + (5 * 15 * p)) / 8;
			constexpr size_t byte_offset_6 = ((m1 - 1) * p + (6 * 15 * p)) / 8;
			constexpr size_t byte_offset_7 = ((m1 - 1) * p + (7 * 15 * p)) / 8;

			// how many bits into the above byte is j + m1*p?
			constexpr size_t bit_offset_0 = ((m1 - 1) * p + (0 * 15 * p)) % 8;
			constexpr size_t bit_offset_1 = ((m1 - 1) * p + (1 * 15 * p)) % 8;
			constexpr size_t bit_offset_2 = ((m1 - 1) * p + (2 * 15 * p)) % 8;
			constexpr size_t bit_offset_3 = ((m1 - 1) * p + (3 * 15 * p)) % 8;
			constexpr size_t bit_offset_4 = ((m1 - 1) * p + (4 * 15 * p)) % 8;
			constexpr size_t bit_offset_5 = ((m1 - 1) * p + (5 * 15 * p)) % 8;
			constexpr size_t bit_offset_6 = ((m1 - 1) * p + (6 * 15 * p)) % 8;
			constexpr size_t bit_offset_7 = ((m1 - 1) * p + (7 * 15 * p)) % 8;

			*(write_t*)(ptr_to_1p + byte_offset_0) &= ~(mask << bit_offset_0);
			*(write_t*)(ptr_to_1p + byte_offset_1) &= ~(mask << bit_offset_1);
			*(write_t*)(ptr_to_1p + byte_offset_2) &= ~(mask << bit_offset_2);
			*(write_t*)(ptr_to_1p + byte_offset_3) &= ~(mask << bit_offset_3);
			*(write_t*)(ptr_to_1p + byte_offset_4) &= ~(mask << bit_offset_4);
			*(write_t*)(ptr_to_1p + byte_offset_5) &= ~(mask << bit_offset_5);
			*(write_t*)(ptr_to_1p + byte_offset_6) &= ~(mask << bit_offset_6);
			*(write_t*)(ptr_to_1p + byte_offset_7) &= ~(mask << bit_offset_7);
		}
		else // two scalar writes required
		{
			clear_next_eight<p, m1>(ptr_to_1p);
			clear_next_eight<p, m2>(ptr_to_1p);
		}
	}



	// for now, assume p % 8 == 1
	void generate_scalar_writes_test_code(size_t p,
										  __attribute__((unused)) size_t j,
										  __attribute__((unused)) uint8_t* const sieve_data)
	{
		//const size_t offset_0 = (j + (1ull * p)) / 8;
		//const size_t offset_1 = (j + (2ull * p)) / 8;
		//const size_t offset_2 = (j + (4ull * p)) / 8;
		//const size_t offset_3 = (j + (7ull * p)) / 8;
		//const size_t offset_4 = (j + (8ull * p)) / 8;
		//const size_t offset_5 = (j + (11ull * p)) / 8;
		//const size_t offset_6 = (j + (13ull * p)) / 8;
		//const size_t offset_7 = (j + (14ull * p)) / 8;

		// prime % 8 yields which of four implementations we need (1, 3, 5, or 7)
		constexpr size_t p_mod_8 = 1;

		std::cout << '\n';

		// generate	64 byte offsets and 64 bit offsets

		std::cout << "Bit offsets for p == " << p << ":\n";
		for (size_t i = 0; i < 8; ++i)
		{
			for (size_t off : { 1, 2, 4, 7, 8, 11, 13, 14 })
			{
				const size_t current = ((i * 15 * p) + ((off - 1) * p)) % 8;

				std::cout << current << ' ';
			}
			std::cout << '\n';
		}
		std::cout << "\n\n\n";

		std::cout << "Bit offsets for p == " << p << ", without using 'p':\n";
		for (size_t i = 0; i < 8; ++i)
		{
			for (size_t off : { 1, 2, 4, 7, 8, 11, 13, 14 })
			{
				const size_t current = ((i * 15 * p_mod_8) + ((off - 1) * p_mod_8)) % 8;

				std::cout << current << ' ';
			}
			std::cout << '\n';
		}
		std::cout << "\n\n\n";

		std::cout << "Byte offsets for p == " << p << ":\n";
		for (size_t i = 0; i < 8; ++i)
		{
			for (size_t off : { 1, 2, 4, 7, 8, 11, 13, 14 })
			{
				const size_t current = ((i * 15 * p) + ((off - 1) * p)) / 8;

				std::cout << current << '\t';
			}
			std::cout << '\n';
		}
		std::cout << "\n\n\n";
	}

}
