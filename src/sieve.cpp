
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



	using sieve_offset_t = util::narrowest_uint_for_val<static_sieve_size>;
	static std::array<sieve_offset_t, small_primes_lookup.size()> sieve_offsets_cache{};

	void set_up_sieve_offsets_cache(const size_t start)
	{
		static_assert(sizeof(sieve_offset_t) >= sizeof(sieve_prime_t));

		// Start with the first prime not in the static sieve.
		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			size_t p = small_primes_lookup[i];

			// We sieve by strides of 15*p, so align p to an (odd) multiple of 15*p
			const size_t p15 = p * 15;

			// Find out how far it is to the next multiple of 15p.
			size_t n = p15 - (start % p15);

			// Start is always odd. Therefore, if n is odd, it is pointing to the next even multiple of p15.
			// -- increase by p15
			if (n % 2 == 1)
				n += p15;
			// However, if n is even, it is pointing to the next odd multiple of p15.
			// -- do nothing

			// handle edge cases where start % prime == 0
			if (n == 2 * p15)
				n = 0;

			// We now have the distance to the next odd multiple of p15.
			// Divide by 2 to get the *index* of the next odd multiple of p15.
			n /= 2;

			if (p <= largest_vector_sieve_prime)
			{
				// We sieve by strides of 8*15*p, starting with a bit offset of 0.
				// Advance by 15*p until we have this alignment.
				while (n % 8 != 0)
					n += p15;

				// If we've ended up at the second multiple of 8*15*p, step back to the first.
				n = util::min(n, n - 8ull * p15);
			}

			sieve_offsets_cache[i] = sieve_offset_t(n);
		}
	}

	void update_sieve_offsets_cache(const uint64_t start,
									const sieve_prime_t* prime_ptr,
									sieve_offset_t* offset_ptr)
	{
		for (size_t prime = *prime_ptr;
			 prime <= largest_sieve_prime;
			 prime = *++prime_ptr, ++offset_ptr)
		{
			// We sieve by strides of 8*15*p, so align p to an (odd) multiple of 8*15*p
			const size_t p15 = prime * 15;

			// Find out how far it is to the next multiple of p15.
			const size_t rem = (start % p15);
			size_t n = p15 - rem;

			// Start is always odd. Therefore:
			// - If n is odd, it is pointing to the next even multiple of p15. Increase by p15.
			// - If n is even, it is pointing to the next odd multiple of p15. Do nothing.
			if (n % 2 == 1)
				n += p15;

			// Handle an edge case where start % prime == 0
			if (rem == 0)
				n = 0;

			// We now have the distance to the next odd multiple of p.
			// Divide by 2 to get the *index* of the next odd multiple of p.
			n /= 2;

			if (prime <= largest_vector_sieve_prime)
			{
				// We sieve by strides of 8*15*p, starting with a bit offset of 0.
				// Advance by 15*p until we have this alignment.
				while (n % 8 != 0)
					n += p15;

				// If we've ended up at the second multiple of 8*15*p, step back to the first.
				n = util::min(n, n - 8ull * p15);
			}

			*offset_ptr = sieve_offset_t(n);
		}
	}

	void detail::verify_sieve_offset_cache(const uint64_t start)
	{
		for (size_t i = static_sieve_primes.size() + 1; small_primes_lookup[i] <= largest_sieve_prime; ++i)
		{
			const uint64_t prime = small_primes_lookup[i];
			const uint64_t p15 = 15ull * prime;
			const uint64_t offset = sieve_offsets_cache[i];

			// 1. start + (2 * offset) should be evenly divisible by 15*p
			if ((start + (2 * offset)) % p15 != 0)
			{
				std::cout << "(start + (2 * offset)) % p15 == " << (start + (2 * offset)) % p15 << ", should be 0\n";
				std::cout << "start == " << start << '\n';
				std::cout << "offset == " << offset << '\n';
				std::cin.ignore();
			}

			// 2. For vectorized sieving, the offset should be evenly divisible by 8.
			if (prime <= largest_vector_sieve_prime && offset % 8 != 0)
			{
				std::cout << "offset % 8 == " << offset % 8 << ", should be 0\n";
				std::cin.ignore();
			}

			// 3. The offset should be less than 8 * 15*p
			if (offset >= 8 * p15)
			{
				std::cout << "prime == " << prime << '\n';
				std::cout << "15*p == " << p15 << '\n';
				std::cout << "8 * 15*p == " << 8 * p15 << '\n';
				std::cout << "offset == " << offset << ", should be < " << 8 * p15 << '\n';
				std::cin.ignore();
			}
		}
	}



	template<size_t p>
	consteval std::array<bit_array<256>, 8> generate_sieve_masks()
	{
		std::array<bit_array<256>, 8> masks{};

		for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset)
		{
			masks[bit_offset].set_all();

			for (size_t i = bit_offset; i < 256; i += p)
			{
				masks[bit_offset].clear_bit(i);
			}
		}

		return masks;
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

	// m_offset = the number of multiples of p past the sieve pointer
	template<size_t p, size_t m_offset = 1>
	__forceinline void generate_vector_writes(uint8_t* const ptr_to_p,
											  const uint256_t& mask_0,
											  const uint256_t& mask_1,
											  const uint256_t& mask_2,
											  const uint256_t& mask_3,
											  const uint256_t& mask_4,
											  const uint256_t& mask_5,
											  const uint256_t& mask_6,
											  const uint256_t& mask_7)
	{
		constexpr size_t stride_size = 8ull * 3 * 5;

		constexpr size_t element_offset = m_offset * p;
		constexpr size_t bit_offset = element_offset % 8;

		constexpr size_t multiples_per_vector_write = 1 + (((256 - 1) - bit_offset) / p);

		if constexpr ((m_offset % 3 == 0 || m_offset % 5 == 0) &&
					  (m_offset + 1) <= stride_size - multiples_per_vector_write + 1)
		{
			// m_offset is divisible by 3 or 5, and we still have room for at least one more write.
			// Increment the offset and continue generating writes.
			generate_vector_writes<p, m_offset + 1>(ptr_to_p, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7);
		}
		else if constexpr (m_offset <= stride_size - multiples_per_vector_write + 1)
		{
			// Generate a vector write for this offset
			constexpr size_t byte_offset = element_offset / 8;

			uint256_t sieve_data = _mm256_loadu_si256((uint256_t*)(ptr_to_p + byte_offset));

			// select from one of eight masks based on the bit offset of this write
			if constexpr (bit_offset == 0) sieve_data = _mm256_and_si256(mask_0, sieve_data);
			if constexpr (bit_offset == 1) sieve_data = _mm256_and_si256(mask_1, sieve_data);
			if constexpr (bit_offset == 2) sieve_data = _mm256_and_si256(mask_2, sieve_data);
			if constexpr (bit_offset == 3) sieve_data = _mm256_and_si256(mask_3, sieve_data);
			if constexpr (bit_offset == 4) sieve_data = _mm256_and_si256(mask_4, sieve_data);
			if constexpr (bit_offset == 5) sieve_data = _mm256_and_si256(mask_5, sieve_data);
			if constexpr (bit_offset == 6) sieve_data = _mm256_and_si256(mask_6, sieve_data);
			if constexpr (bit_offset == 7) sieve_data = _mm256_and_si256(mask_7, sieve_data);

			_mm256_storeu_si256((uint256_t*)(ptr_to_p + byte_offset), sieve_data);

			// Advance by the size of the write and continue generating writes.
			generate_vector_writes<p, m_offset + multiples_per_vector_write>(ptr_to_p, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7);
		}
	}

	template<size_t p>
	__forceinline void vectorized_sieve_pass(sieve_container& sieve,
											 const sieve_prime_t*& prime_ptr,
											 sieve_offset_t*& offset_cache_ptr)
	{
		constexpr size_t sieve_end = sieve_container::size();
		constexpr size_t padded_end = sieve_end - (8 * 15 * p);

		constexpr static std::array<bit_array<256>, 8> sieve_masks = generate_sieve_masks<p>();

		// for every 15 multiples of a prime p:
		//
		// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
		//    x  x     x        x  x        x     x  x   <-- values to mark composite/false
		// x        x     x  x        x  x     x         <-- values to ignore

		// Get the position of the next odd multiple of 8*15*p
		size_t j = *offset_cache_ptr;

		const uint256_t mask_0 = _mm256_loadu_si256((const uint256_t*)sieve_masks[0].data());
		const uint256_t mask_1 = _mm256_loadu_si256((const uint256_t*)sieve_masks[1].data());
		const uint256_t mask_2 = _mm256_loadu_si256((const uint256_t*)sieve_masks[2].data());
		const uint256_t mask_3 = _mm256_loadu_si256((const uint256_t*)sieve_masks[3].data());
		const uint256_t mask_4 = _mm256_loadu_si256((const uint256_t*)sieve_masks[4].data());
		const uint256_t mask_5 = _mm256_loadu_si256((const uint256_t*)sieve_masks[5].data());
		const uint256_t mask_6 = _mm256_loadu_si256((const uint256_t*)sieve_masks[6].data());
		const uint256_t mask_7 = _mm256_loadu_si256((const uint256_t*)sieve_masks[7].data());

		{
			// We've made sure j is evenly divisible by 8
			uint8_t* sieve_ptr = sieve.data() + (j / 8);

			do
			{
				// 89 is the first prime that would require at least one final scalar write.
				// We don't handle this right now.
				static_assert(p < 89);

				// Sieve by big strides of 8*15 * p.
				// This is guaranteed to have a period of one (ie, the first bit offset is always the same)
				generate_vector_writes<p>(sieve_ptr, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7);

				j += 8 * 15 * p;
				sieve_ptr += (8 * 15 * p) / 8;

			} while (j <= padded_end);
		}

		// If we run the hot loop with a bit offset of 0, the size of the sieve will
		// introduce the same new offset per iteration.
		// For a given offset_increase, how many multiples of 15*p is it to the next bit offset of 0?
		constexpr size_t step = []() consteval {
			constexpr size_t new_offset_per_iter = 8ull - (sieve_container::size() % 8);
			size_t n = new_offset_per_iter;
			size_t multiples_of_15p = 0;
			while (n % 8 != 0)
			{
				n += 15 * p;
				++multiples_of_15p;
			}

			return multiples_of_15p;
		}();

		// Calculate and cache the offset for the next sieving.

		// += 8 * 15*p advances one full stride, taking j past the end of the sieve. j % 8 is still 0.
		// += step * 15*p cancels out the %8 offset introduced when we subtract the size of the sieve.
		j += (8ull + step) * 15 * p;

		// The distance past the end of the sieve is the index for the start of the next sieving.
		j -= sieve_container::size();

		// If we've ended up at the second multiple of 8*15*p, step back to the first.
		j = util::min(j, j - 8ull * 15 * p);

		*offset_cache_ptr = sieve_offset_t(j);

		++prime_ptr;
		++offset_cache_ptr;
	}

	__forceinline size_t clear_bit(uint8_t* const sieve_data,
								   size_t offset, size_t delta_to_next_offset)
	{
		size_t next_offset = offset + delta_to_next_offset;

		size_t bit = _pdep_u64(offset, 7);
		offset >>= 3;
		sieve_data[offset] &= ~(1 << bit);

		return next_offset;
	}

	__forceinline size_t clear_8_of_15(const size_t p, size_t j, uint8_t* const sieve_data)
	{
		j += p;
		j = clear_bit(sieve_data, j, 1 * p); // mark j + p, advance by 1p
		j = clear_bit(sieve_data, j, 2 * p); // mark j + 2p, advance by 2p
		j = clear_bit(sieve_data, j, 3 * p); // mark j + 4p, advance by 3p
		j = clear_bit(sieve_data, j, 1 * p); // mark j + 7p, advance by 1p
		j = clear_bit(sieve_data, j, 3 * p); // mark j + 8p, advance by 3p
		j = clear_bit(sieve_data, j, 2 * p); // mark j + 11p, advance by 2p
		j = clear_bit(sieve_data, j, 1 * p); // mark j + 13p, advance by 1p
		j = clear_bit(sieve_data, j, 1 * p); // mark j + 14p, advance by 1p

		return j;
	}

	__forceinline size_t clear_bit(uint8_t* const sieve_data,
								   size_t offset, size_t delta_to_next_offset,
								   const size_t mask)
	{
		size_t next_offset = offset + delta_to_next_offset;

		offset >>= 3;
		sieve_data[offset] &= mask;

		return next_offset;
	}

	__forceinline size_t clear_8_of_15(const size_t p, size_t j, uint8_t* const sieve_data,
									   const size_t mask_a,
									   const size_t mask_b,
									   const size_t mask_c,
									   const size_t mask_d,
									   const size_t mask_e,
									   const size_t mask_f,
									   const size_t mask_g,
									   const size_t mask_h)
	{
		j += p;
		j = clear_bit(sieve_data, j, 1 * p, mask_a); // mark j + p, advance by 1p
		j = clear_bit(sieve_data, j, 2 * p, mask_b); // mark j + 2p, advance by 2p
		j = clear_bit(sieve_data, j, 3 * p, mask_c); // mark j + 4p, advance by 3p
		j = clear_bit(sieve_data, j, 1 * p, mask_d); // mark j + 7p, advance by 1p
		j = clear_bit(sieve_data, j, 3 * p, mask_e); // mark j + 8p, advance by 3p
		j = clear_bit(sieve_data, j, 2 * p, mask_f); // mark j + 11p, advance by 2p
		j = clear_bit(sieve_data, j, 1 * p, mask_g); // mark j + 13p, advance by 1p
		j = clear_bit(sieve_data, j, 1 * p, mask_h); // mark j + 14p, advance by 1p

		return j;
	}

	void partial_sieve(const uint64_t number,
					   sieve_container& sieve,
					   const size_t sieve_popcount)
	{
		// Sieve primes by strides of 15*p:
		//
		// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
		//    x  x     x        x  x        x     x  x   <-- values to mark composite/false
		// x        x     x  x        x  x     x         <-- values to ignore

		constexpr size_t sieve_end = sieve_container::size();

		// Start with the first prime not in the static sieve
		const sieve_prime_t* prime_ptr = small_primes_lookup.data() + static_sieve_primes.size() + 1;
		sieve_offset_t* offset_cache_ptr = sieve_offsets_cache.data() + static_sieve_primes.size() + 1;

		// calculate sieve density
		double density = double(sieve_popcount) / sieve_container::size();

		// don't do any sieving if our bit pattern filters + static sieve already cleared enough
		if (density < vector_density_threshold)
		{
			update_sieve_offsets_cache(number + 2 * sieve_container::size(), prime_ptr, offset_cache_ptr);
			return;
		}

		static_assert(small_primes_lookup[static_sieve_primes.size() + 1] == 19);
		vectorized_sieve_pass<19>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<23>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<29>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<31>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<37>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<41>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<43>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<47>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<53>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<59>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<61>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<67>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<71>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<73>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<79>(sieve, prime_ptr, offset_cache_ptr);
		static_assert(largest_vector_sieve_prime == 79);

		constexpr double vectorized_sieving_removes = .3106; // 31.06%
		constexpr double scale = 1.0 - vectorized_sieving_removes;

		density *= scale;

		for (;;)
		{
			if (density < scalar_density_threshold)
			{
				// If we stop sieving early, we still need to update our offsets cache
				update_sieve_offsets_cache(number + 2 * sieve_container::size(),
										   prime_ptr, offset_cache_ptr);
				break;
			}

			// Get the next prime
			const size_t p = *prime_ptr++;

			// Update the density estimate
			density -= (1.0 / p) * density;

			// Get the position of the next odd multiple of p*15
			size_t j = *offset_cache_ptr;

			// Stop marking 15*p early (don't handle padding)
			const size_t padded_end = sieve_end - (15 * p);

			const size_t extra_padded_end = sieve_end - (15ull * 8 * p);

			uint8_t* const sieve_data = sieve.data();

			// run 0-7 steps to align the hot loop's first write to a bit offset of 0
			while ((j + p) % 8 != 0)
			{
				j = clear_8_of_15(p, j, sieve_data);
			}

			// generate eight byte masks in usage order, starting from (j + p) % 8 == 0
			const size_t mask_a = ~(1ull << ((j + (1ull * p)) % 8));
			const size_t mask_b = ~(1ull << ((j + (2ull * p)) % 8));
			const size_t mask_c = ~(1ull << ((j + (4ull * p)) % 8));
			const size_t mask_d = ~(1ull << ((j + (7ull * p)) % 8));
			const size_t mask_e = ~(1ull << ((j + (8ull * p)) % 8));
			const size_t mask_f = ~(1ull << ((j + (11ull * p)) % 8));
			const size_t mask_g = ~(1ull << ((j + (13ull * p)) % 8));
			const size_t mask_h = ~(1ull << ((j + (14ull * p)) % 8));

			// Each iteration, j points to an (already marked) multiple of p*15.
			// Clear eight multiples of p within each of eight strides of p*15.
			// The first write in the first stride has a bit offset of 0.
			do
			{
				j = clear_8_of_15(p, j, sieve_data, mask_a, mask_b, mask_c, mask_d, mask_e, mask_f, mask_g, mask_h);
				j = clear_8_of_15(p, j, sieve_data, mask_e, mask_a, mask_f, mask_h, mask_d, mask_b, mask_c, mask_g);
				j = clear_8_of_15(p, j, sieve_data, mask_d, mask_e, mask_b, mask_g, mask_h, mask_a, mask_f, mask_c);
				j = clear_8_of_15(p, j, sieve_data, mask_h, mask_d, mask_a, mask_c, mask_g, mask_e, mask_b, mask_f);
				j = clear_8_of_15(p, j, sieve_data, mask_g, mask_h, mask_e, mask_f, mask_c, mask_d, mask_a, mask_b);
				j = clear_8_of_15(p, j, sieve_data, mask_c, mask_g, mask_d, mask_b, mask_f, mask_h, mask_e, mask_a);
				j = clear_8_of_15(p, j, sieve_data, mask_f, mask_c, mask_h, mask_a, mask_b, mask_g, mask_d, mask_e);
				j = clear_8_of_15(p, j, sieve_data, mask_b, mask_f, mask_g, mask_e, mask_a, mask_c, mask_h, mask_d);
			} while (j < extra_padded_end);

			// cleanup loop
			while (j < padded_end)
			{
				j = clear_8_of_15(p, j, sieve_data);
			}

			// Calculate and cache the offset for the next sieving
			*offset_cache_ptr++ = sieve_offset_t((j + (15 * p)) - sieve_end);

			if (p == largest_sieve_prime) break;
		}
	}



	namespace detail
	{
		constexpr uint64_t n_sieve_chunks = (sieve_container::size() / 64) + (sieve_container::size() % 64 != 0);
		constexpr size_t max_pc = 27; // 27 is the highest popcount of 64-bit chunks in the static sieve

		using chunk_count_t = util::narrowest_uint_for_val<n_sieve_chunks>;
		using chunk_idx_t = chunk_count_t;

		// [popcount] -> number of chunks with popcount
		static std::array<chunk_count_t, max_pc + 1> n_chunks_with_pc;
		// [popcount][chunks]
		static std::array<std::array<uint64_t, n_sieve_chunks>, max_pc + 1> sorted_chunks;
		// [popcount][chunk indexes]
		static std::array<std::array<chunk_count_t, n_sieve_chunks>, max_pc + 1> chunk_indexes;

		template<size_t popcount>
		__forceinline void extract_candidates(uint64_t& chunk,
											  const uint64_t chunk_index,
											  uint64_t*& candidates)
		{
			// read the next bit in the chunk
			const size_t tzcnt = _tzcnt_u64(chunk);

			// reset the bit we just read
			if constexpr (popcount > 1)
				chunk = _blsr_u64(chunk);

			// store the chunk's index plus the bit's index
			*candidates++ = chunk_index + tzcnt;

			if constexpr (popcount - 1 > 0)
				extract_candidates<popcount - 1>(chunk, chunk_index, candidates);
		}

		template<size_t popcount>
		__forceinline void extract_candidates_with_popcount(uint64_t*& candidates)
		{
			for (size_t j = 0, n_chunks = n_chunks_with_pc[popcount]; j < n_chunks; ++j)
			{
				uint64_t chunk = sorted_chunks[popcount][j];
				uint64_t index = chunk_indexes[popcount][j] * 64ull;

				// generate instructions to extract n candidates
				extract_candidates<popcount>(chunk, index, candidates);
			}
		}

		__forceinline void sort_chunks(const uint64_t* sieve_data)
		{
			using namespace detail;

			constexpr size_t n_chunks_rounded = n_sieve_chunks - (n_sieve_chunks % 3);

			// load one iteration ahead
			uint64_t chunk_0 = sieve_data[0];
			uint64_t chunk_1 = sieve_data[1];
			uint64_t chunk_2 = sieve_data[2];
			for (size_t i = 0; i < n_chunks_rounded; i += 3)
			{
				const size_t pc_0 = pop_count(chunk_0);
				const size_t pc_1 = pop_count(chunk_1);
				const size_t pc_2 = pop_count(chunk_2);

				const chunk_count_t idx_0 = n_chunks_with_pc[pc_0]++;
				const chunk_count_t idx_1 = n_chunks_with_pc[pc_1]++;
				const chunk_count_t idx_2 = n_chunks_with_pc[pc_2]++;

				sorted_chunks[pc_0][idx_0] = chunk_0;
				chunk_indexes[pc_0][idx_0] = chunk_idx_t(i + 0);
				sorted_chunks[pc_1][idx_1] = chunk_1;
				chunk_indexes[pc_1][idx_1] = chunk_idx_t(i + 1);
				sorted_chunks[pc_2][idx_2] = chunk_2;
				chunk_indexes[pc_2][idx_2] = chunk_idx_t(i + 2);

				chunk_0 = sieve_data[i + 3 + 0];
				chunk_1 = sieve_data[i + 3 + 1];
				chunk_2 = sieve_data[i + 3 + 2];
			}

			for (size_t i = n_chunks_rounded; i < n_sieve_chunks; ++i)
			{
				const uint64_t chunk = sieve_data[i];

				const size_t pc = pop_count(chunk);
				const size_t idx = n_chunks_with_pc[pc]++;

				sorted_chunks[pc][idx] = chunk;
				chunk_indexes[pc][idx] = chunk_idx_t(i);
			}
		}

		__forceinline void sort_chunks_vectorized(const uint64_t* sieve_data)
		{
			using namespace detail;

			constexpr static uint256_t static_identity = { .m256i_u32{ 0, 1, 2, 3, 4, 5, 6, 7 } };
			constexpr static uint256_t static_inc = { .m256i_u32{ 1, 1, 1, 1, 1, 1, 1, 1 } };

			const uint256_t identity = _mm256_loadu_si256(&static_identity);
			const uint256_t inc = _mm256_loadu_si256(&static_inc);
			const uint256_t zero_register = _mm256_setzero_si256();

			uint256_t pc_counts{};

			for (size_t i = 0; i < n_sieve_chunks; ++i)
			{
				const uint64_t chunk = sieve_data[i];
				const uint32_t pc = (uint32_t)pop_count(chunk);
				uint32_t idx = 0;

				if (pc < 8)
				{
					const uint256_t bcast_pc = _mm256_set1_epi32(pc); // copy the popcount (index) to all positions
					const uint256_t mask = _mm256_cmpeq_epi32(bcast_pc, identity); // mask to select the element at [pc]
					const uint256_t masked_inc = _mm256_and_si256(mask, inc); // select one increment

					const uint256_t shuffled = _mm256_permutevar8x32_epi32(pc_counts, bcast_pc); // get existing idx

					pc_counts = _mm256_add_epi16(pc_counts, masked_inc); // add one to the pc count

					idx = _mm256_extract_epi32(shuffled, 0); // extract idx
				}
				else
				{
					idx = n_chunks_with_pc[pc]++;
				}

				sorted_chunks[pc][idx] = chunk;
				chunk_indexes[pc][idx] = chunk_idx_t(i);
			}

			for (size_t i = 0; i < 8; ++i)
			{
				n_chunks_with_pc[i] = chunk_count_t(pc_counts.m256i_u32[i]);
			}
		}

		__forceinline void extract_candidates_as_offsets(uint64_t*& candidates)
		{
			using namespace detail;

			extract_candidates_with_popcount<1>(candidates);
			extract_candidates_with_popcount<2>(candidates);
			extract_candidates_with_popcount<3>(candidates);
			extract_candidates_with_popcount<4>(candidates);
			extract_candidates_with_popcount<5>(candidates);
			extract_candidates_with_popcount<6>(candidates);
			extract_candidates_with_popcount<7>(candidates);
			extract_candidates_with_popcount<8>(candidates);
			extract_candidates_with_popcount<9>(candidates);
			extract_candidates_with_popcount<10>(candidates);
			extract_candidates_with_popcount<11>(candidates);
			extract_candidates_with_popcount<12>(candidates);
			extract_candidates_with_popcount<13>(candidates);
			extract_candidates_with_popcount<14>(candidates);
			extract_candidates_with_popcount<15>(candidates);
			extract_candidates_with_popcount<16>(candidates);
			extract_candidates_with_popcount<17>(candidates);
			extract_candidates_with_popcount<18>(candidates);
			extract_candidates_with_popcount<19>(candidates);
			extract_candidates_with_popcount<20>(candidates);
			extract_candidates_with_popcount<21>(candidates);
			extract_candidates_with_popcount<22>(candidates);
			extract_candidates_with_popcount<23>(candidates);
			extract_candidates_with_popcount<24>(candidates);
			extract_candidates_with_popcount<25>(candidates);
			extract_candidates_with_popcount<26>(candidates);
			extract_candidates_with_popcount<27>(candidates);
			static_assert(max_pc == 27);
		}

		__forceinline void convert_candidates_to_bitstrings(uint64_t* const candidates_start,
															const uint64_t* const candidates_end,
															const uint64_t number)
		{
			// autovectorizes
			for (uint64_t* ptr = candidates_start; ptr < candidates_end; ++ptr)
			{
				uint64_t candidate = *ptr;
				candidate *= 2;
				candidate += number;
				*ptr = candidate;
			}
		}
	}

	// sort 64-bit chunks of the sieve by popcount, then use custom loops that perform the exact number of required reads
	uint64_t* gather_sieve_results(uint64_t* candidates,
								   const sieve_container& sieve,
								   const uint64_t number)
	{
		using namespace detail;

		// autovectorizes to three large writes
		for (chunk_count_t& pc : n_chunks_with_pc)
			pc = 0;

		// prepare arrays of 64-bit chunks, where the chunks in array n contain exactly n set bits (n candidates)
		sort_chunks((const uint64_t*)sieve.data());

		uint64_t* const candidates_start = candidates;

		extract_candidates_as_offsets(candidates);

		convert_candidates_to_bitstrings(candidates_start, candidates, number);

		return candidates;
	}

}
