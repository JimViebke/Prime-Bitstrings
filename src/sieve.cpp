
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
			sieve_offset_t p = sieve_offset_t(small_primes_lookup[i]);

			// We sieve by strides of 15*p, so align p to an (odd) multiple of 15*p
			p *= 15;

			// Find out how far it is to the next multiple of p.
			sieve_offset_t n = p - (start % p);

			// Start is always odd. Therefore, if n is odd, it is pointing to the next even multiple of p.
			// -- increase by p
			if (n % 2 == 1)
				n += p;
			// However, if n is even, it is pointing to the next odd multiple of p.
			// -- do nothing

			// handle edge cases where start % prime == 0
			if (n == 2 * p)
				n = 0;

			// We now have the distance to the next odd multiple of p.
			// Divide by 2 to store the *index* of the next odd multiple of p.
			sieve_offsets_cache[i] = n / 2;
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
			// We sieve by strides of 15*p, so align p to an (odd) multiple of 15*p
			prime *= 15;

			// Find out how far it is to the next multiple of p.
			const size_t rem = (start % prime);
			size_t n = prime - rem;

			// Start is always odd. Therefore:
			// - If n is odd, it is pointing to the next even multiple of p. Increase by p.
			// - If n is even, it is pointing to the next odd multiple of p. Do nothing.
			if (n % 2 == 1)
				n += prime;

			// Handle an edge case where start % prime == 0
			if (rem == 0)
				n = 0;

			// Divide by 2 to convert the distance to the next odd multiple of p
			// to the *index* of the next odd multiple of p.
			*offset_ptr = sieve_offset_t(n / 2);
		}
	}

	void detail::verify_sieve_offset_cache(const uint64_t start)
	{
		// - start + (2 * offset) should be evenly divisible by 15*p
		// - the offset should be less than 15*p

		for (size_t i = static_sieve_primes.size() + 1; small_primes_lookup[i] <= largest_sieve_prime; ++i)
		{
			const uint64_t p15 = 15ull * small_primes_lookup[i];
			const uint64_t j = sieve_offsets_cache[i];

			const uint64_t remainder = (start + (2 * j)) % p15;

			if (j >= p15 || remainder != 0)
			{
				std::cout << "small_primes_lookup[" << i << "] = " << small_primes_lookup[i] << '\n';
				std::cout << "sieve_offsets_cache[" << i << "] = " << j << '\n';
				std::cout << "max allowed offset (p*15) = " << p15 << '\n';
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

	// for a write to j + 4*p: eight_vector_writes<p, 4-1>(ptr_to_1p)
	template<size_t p, size_t m_offset>
	__forceinline void eight_vector_writes(uint8_t* ptr_to_1p,
										   const uint256_t& mask_0,
										   const uint256_t& mask_1,
										   const uint256_t& mask_2,
										   const uint256_t& mask_3,
										   const uint256_t& mask_4,
										   const uint256_t& mask_5,
										   const uint256_t& mask_6,
										   const uint256_t& mask_7)
	{
		constexpr size_t byte_index_0 = ((m_offset * p) + (0 * 15 * p)) / 8; // always 0 when m_offset is 0
		constexpr size_t byte_index_1 = ((m_offset * p) + (1 * 15 * p)) / 8;
		constexpr size_t byte_index_2 = ((m_offset * p) + (2 * 15 * p)) / 8;
		constexpr size_t byte_index_3 = ((m_offset * p) + (3 * 15 * p)) / 8;
		constexpr size_t byte_index_4 = ((m_offset * p) + (4 * 15 * p)) / 8;
		constexpr size_t byte_index_5 = ((m_offset * p) + (5 * 15 * p)) / 8;
		constexpr size_t byte_index_6 = ((m_offset * p) + (6 * 15 * p)) / 8;
		constexpr size_t byte_index_7 = ((m_offset * p) + (7 * 15 * p)) / 8;

		uint256_t sieve_data_0 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_0));
		uint256_t sieve_data_1 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_1));
		uint256_t sieve_data_2 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_2));
		uint256_t sieve_data_3 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_3));
		uint256_t sieve_data_4 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_4));
		uint256_t sieve_data_5 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_5));
		uint256_t sieve_data_6 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_6));
		uint256_t sieve_data_7 = _mm256_loadu_si256((uint256_t*)(ptr_to_1p + byte_index_7));

		sieve_data_0 = _mm256_and_si256(mask_0, sieve_data_0);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_0), sieve_data_0);

		sieve_data_1 = _mm256_and_si256(mask_1, sieve_data_1);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_1), sieve_data_1);

		sieve_data_2 = _mm256_and_si256(mask_2, sieve_data_2);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_2), sieve_data_2);

		sieve_data_3 = _mm256_and_si256(mask_3, sieve_data_3);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_3), sieve_data_3);

		sieve_data_4 = _mm256_and_si256(mask_4, sieve_data_4);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_4), sieve_data_4);

		sieve_data_5 = _mm256_and_si256(mask_5, sieve_data_5);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_5), sieve_data_5);

		sieve_data_6 = _mm256_and_si256(mask_6, sieve_data_6);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_6), sieve_data_6);

		sieve_data_7 = _mm256_and_si256(mask_7, sieve_data_7);
		_mm256_storeu_si256((uint256_t*)(ptr_to_1p + byte_index_7), sieve_data_7);
	}

	template<size_t p>
	consteval size_t hi_mask_offset()
	{
		// begin at the byte containing 8p
		if constexpr (p == 37 || p == 41 || p == 43 || p == 47)
			return 8;

		// begin at the byte containing 7p
		if constexpr (p == 53 || p == 59 || p == 61)
			return 7;

		// begin at the byte containing 11p
		if constexpr (p == 67 || p == 71 || p == 73 || p == 79)
			return 11;

		// compile-time error if we don't find an answer
	}

	template<size_t p>
	__forceinline void vectorized_sieve_pass(sieve_container& sieve,
											 const sieve_prime_t*& prime_ptr,
											 sieve_offset_t*& offset_cache_ptr)
	{
		constexpr size_t sieve_end = sieve_container::size();
		constexpr size_t padded_end = sieve_end - (15 * p);

		constexpr static std::array<bit_array<256>, 8> sieve_masks = generate_sieve_masks<p>();

		// Get the position of the next odd multiple of p*15
		size_t j = *offset_cache_ptr;

		// 0-7 leading steps to get the hot loop to start with bit_idx == 0
		for (size_t bit_index = 0; (bit_index = (j + p) % 8) != 0; )
		{
			uint256_t* const sieve_ptr = (uint256_t*)(sieve.data() + (j + p) / 8);

			const uint256_t mask = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index].data());
			uint256_t sieve_data = _mm256_loadu_si256(sieve_ptr);

			if constexpr (p >= 37)
			{
				const size_t byte_index_hi = (j + (hi_mask_offset<p>() * p)) / 8;
				const size_t bit_index_hi = (j + (hi_mask_offset<p>() * p)) % 8;

				uint256_t* const sieve_ptr_hi = (uint256_t*)(sieve.data() + byte_index_hi);

				const uint256_t mask_hi = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index_hi].data());
				uint256_t sieve_data_hi = _mm256_loadu_si256(sieve_ptr_hi);

				sieve_data = _mm256_and_si256(mask, sieve_data);
				sieve_data_hi = _mm256_and_si256(mask_hi, sieve_data_hi);
				_mm256_storeu_si256(sieve_ptr, sieve_data);
				_mm256_storeu_si256(sieve_ptr_hi, sieve_data_hi);
			}
			else // everyone else
			{
				sieve_data = _mm256_and_si256(mask, sieve_data);
				_mm256_storeu_si256(sieve_ptr, sieve_data);
			}

			if constexpr (p == 29 || p == 31)
			{
				sieve.clear_bit(j + 11 * p);
			}

			if constexpr (p == 23 || p == 29 || p == 31 || p == 53 || p == 59 || p == 61)
			{
				clear_two_adjacent<p>(sieve, j + 13 * p);
			}

			if constexpr (p == 43 || p == 47)
			{
				sieve.clear_bit(j + (7 * p));
				sieve.clear_bit(j + (14 * p));
			}

			if constexpr (p == 67 || p == 71 || p == 73 || p == 79)
			{
				clear_two_adjacent<p>(sieve, j + 7 * p);
			}

			j += (15 * p);
		}

		constexpr size_t extra_padded_end = sieve_end - (8 * 15 * p);

		// We've made sure j+p is evenly divisible by 8, so we can remove j+p from
		// our offset calculations
		const uint256_t mask_0 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(0 * 15 * p) % 8].data());
		const uint256_t mask_1 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(1 * 15 * p) % 8].data());
		const uint256_t mask_2 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(2 * 15 * p) % 8].data());
		const uint256_t mask_3 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(3 * 15 * p) % 8].data());
		const uint256_t mask_4 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(4 * 15 * p) % 8].data());
		const uint256_t mask_5 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(5 * 15 * p) % 8].data());
		const uint256_t mask_6 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(6 * 15 * p) % 8].data());
		const uint256_t mask_7 = _mm256_loadu_si256((const uint256_t*)sieve_masks[(7 * 15 * p) % 8].data());

		{
			uint8_t* sieve_ptr = sieve.data() + ((j + p) / 8);

			do
			{
				eight_vector_writes<p, 0>(sieve_ptr, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7);

				if constexpr (p >= 37)
				{
					// We have eight possible orders that the high masks can be read in.
					// At compile time, check which of eight masks has the required bit_offset,
					//  then perform the eight high writes starting with that mask.
					constexpr size_t p_offset = hi_mask_offset<p>() - 1;
					constexpr size_t high_bit_offset = (p_offset * p) % 8;

					if constexpr (high_bit_offset == (0 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7);
					if constexpr (high_bit_offset == (1 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7, mask_0);
					if constexpr (high_bit_offset == (2 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_2, mask_3, mask_4, mask_5, mask_6, mask_7, mask_0, mask_1);
					if constexpr (high_bit_offset == (3 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_3, mask_4, mask_5, mask_6, mask_7, mask_0, mask_1, mask_2);
					if constexpr (high_bit_offset == (4 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_4, mask_5, mask_6, mask_7, mask_0, mask_1, mask_2, mask_3);
					if constexpr (high_bit_offset == (5 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_5, mask_6, mask_7, mask_0, mask_1, mask_2, mask_3, mask_4);
					if constexpr (high_bit_offset == (6 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_6, mask_7, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5);
					if constexpr (high_bit_offset == (7 * 15 * p) % 8)
						eight_vector_writes<p, p_offset>(sieve_ptr, mask_7, mask_0, mask_1, mask_2, mask_3, mask_4, mask_5, mask_6);
				}

				if constexpr (p == 29 || p == 31)
				{
					clear_next_eight<p, 11>(sieve_ptr);
				}

				if constexpr (p == 23 || p == 29 || p == 31 || p == 53 || p == 59 || p == 61)
				{
					clear_next_eight<p, 13, 14>(sieve_ptr);
				}

				if constexpr (p == 43 || p == 47)
				{
					clear_next_eight<p, 7>(sieve_ptr);
					clear_next_eight<p, 14>(sieve_ptr);
				}

				if constexpr (p == 67 || p == 71 || p == 73 || p == 79)
				{
					clear_next_eight<p, 7, 8>(sieve_ptr);
				}

				j += 8 * 15 * p;
				sieve_ptr += (8 * 15 * p) / 8;

			} while (j < extra_padded_end);
		}

		// cleanup loop
		while (j < padded_end) // stop marking 15*p early (don't handle padding)
		{
			const size_t bit_index = (j + p) % 8;

			uint256_t* const sieve_ptr = (uint256_t*)(sieve.data() + (j + p) / 8);

			const uint256_t mask = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index].data());
			uint256_t sieve_data = _mm256_loadu_si256(sieve_ptr);

			if constexpr (p >= 37)
			{
				const size_t byte_index_hi = (j + (hi_mask_offset<p>() * p)) / 8;
				const size_t bit_index_hi = (j + (hi_mask_offset<p>() * p)) % 8;

				uint256_t* const sieve_ptr_hi = (uint256_t*)(sieve.data() + byte_index_hi);

				const uint256_t mask_hi = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index_hi].data());
				uint256_t sieve_data_hi = _mm256_loadu_si256(sieve_ptr_hi);

				sieve_data = _mm256_and_si256(mask, sieve_data);
				sieve_data_hi = _mm256_and_si256(mask_hi, sieve_data_hi);
				_mm256_storeu_si256(sieve_ptr, sieve_data);
				_mm256_storeu_si256(sieve_ptr_hi, sieve_data_hi);
			}
			else // everyone else
			{
				sieve_data = _mm256_and_si256(mask, sieve_data);
				_mm256_storeu_si256(sieve_ptr, sieve_data);
			}

			if constexpr (p == 29 || p == 31)
			{
				sieve.clear_bit(j + 11 * p);
			}

			if constexpr (p == 23 || p == 29 || p == 31 || p == 53 || p == 59 || p == 61)
			{
				clear_two_adjacent<p>(sieve, j + 13 * p);
			}

			if constexpr (p == 43 || p == 47)
			{
				sieve.clear_bit(j + (7 * p));
				sieve.clear_bit(j + (14 * p));
			}

			if constexpr (p == 67 || p == 71 || p == 73 || p == 79)
			{
				clear_two_adjacent<p>(sieve, j + 7 * p);
			}

			j += (15 * p);
		}

		// Calculate and cache the offset for the next sieving
		*offset_cache_ptr = sieve_offset_t((j + (15 * p)) - sieve_end);

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
					   sieve_container& sieve)
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
		double density = double(sieve.count_bits()) / sieve_container::size();

		// don't do any sieving if our bit pattern filters + static sieve already cleared enough
		if (density < vector_density_threshold)
		{
			update_sieve_offsets_cache(number + 2 * sieve_container::size(), prime_ptr, offset_cache_ptr);
			return;
		}

		if constexpr (small_primes_lookup[static_sieve_primes.size() + 1] == 19)
		{
			// one simd write from 1p, plus 0-1 scalar writes for 11, and 0-1 scalar writes for 13,14
			vectorized_sieve_pass<19>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<23>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<29>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<31>(sieve, prime_ptr, offset_cache_ptr);

			// two simd writes from 1p and 8p, plus 0-2 scalar writes for 7,14
			vectorized_sieve_pass<37>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<41>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<43>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<47>(sieve, prime_ptr, offset_cache_ptr);

			// two simd writes from 1p and 7p, plus 0-2 scalar writes for 13,14
			vectorized_sieve_pass<53>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<59>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<61>(sieve, prime_ptr, offset_cache_ptr);

			// two simd writes from 1p and 11p, plus 2 scalar writes for 7,8
			vectorized_sieve_pass<67>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<71>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<73>(sieve, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<79>(sieve, prime_ptr, offset_cache_ptr);
		}

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
		static std::array<chunk_count_t, max_pc + 1> n_chunks_with_pc;
		// [popcount][chunk]
		static std::array<std::array<uint64_t, n_sieve_chunks>, max_pc + 1> sorted_chunks;

		using chunk_idx_t = util::narrowest_uint_for_val<n_sieve_chunks>;
		// [popcount][idx]
		static std::array<std::array<chunk_idx_t, n_sieve_chunks>, max_pc + 1> chunk_indexes;

		template<size_t n_bits>
		__forceinline void extract_candidates(uint64_t& chunk,
											  const uint64_t idx,
											  uint64_t*& candidates)
		{
			// read the next bit in the chunk
			const size_t tzcnt = _tzcnt_u64(chunk);

			// reset the bit we just read
			if constexpr (n_bits > 1)
				chunk = _blsr_u64(chunk);

			// store the chunk's index plus the bit's index
			*candidates++ = idx + tzcnt;

			if constexpr (n_bits - 1 > 0)
				extract_candidates<n_bits - 1>(chunk, idx, candidates);
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
	}

	// sort 64-bit chunks of the sieve by popcount, then use custom loops that perform the exact number of required reads
	uint64_t* gather_sieve_results(uint64_t* candidates,
								   const sieve_container& sieve,
								   const uint64_t number)
	{
		using namespace detail;

		const uint64_t* sieve_data = (uint64_t*)sieve.data();

		// autovectorizes to three large writes
		for (chunk_count_t& pc : n_chunks_with_pc)
			pc = 0;

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

			const size_t idx_0 = n_chunks_with_pc[pc_0]++;
			const size_t idx_1 = n_chunks_with_pc[pc_1]++;
			const size_t idx_2 = n_chunks_with_pc[pc_2]++;

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

		uint64_t* const candidates_start = candidates;

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

		// autovectorizes
		for (uint64_t* ptr = candidates_start; ptr < candidates; ++ptr)
		{
			uint64_t candidate = *ptr;
			candidate *= 2;
			candidate += number;
			*ptr = candidate;
		}

		return candidates;
	}

}
