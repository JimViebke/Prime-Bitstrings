#pragma once

#include "util/types.hpp"

namespace mbp
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

	// Only step by 15 using primes that make at least X writes per sieving, where X is a tuneable knob.
	consteval sieve_prime_t get_threshold(const size_t X)
	{
		for (size_t i = 1; i < small_primes_lookup.size(); ++i)
			if (small_primes_lookup[i] > (static_sieve_size / X))
				return small_primes_lookup[i - 1];
		// compile-time error if we don't find a valid answer
	}
	constexpr sieve_prime_t last_prime_for_stepping_by_fifteen = 263; // get_threshold((5 + 1) * 15);

	using sieve_offset_t = util::narrowest_uint_for_val<static_sieve_size>;
	std::vector<sieve_offset_t> sieve_offsets_cache(small_primes_lookup.size());

	void set_up_sieve_offsets_cache(const size_t start)
	{
		static_assert(sizeof(sieve_offset_t) >= sizeof(sieve_prime_t));

		// Start with the first prime not in the static sieve.
		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			sieve_offset_t p = sieve_offset_t(small_primes_lookup[i]);

			// p has stricter "alignment" requirements when sieving in steps of p*n
			if (p <= last_prime_for_stepping_by_fifteen)
			{
				p *= 15;
			}
			else // step by threes
			{
				p *= 3;
			}

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
			 prime <= last_prime_for_stepping_by_fifteen;
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

	void verify_sieve_offset_cache(const uint64_t start)
	{
		// - start + (2 * offset) should be evenly divisible by 15*p
		// - the offset should be less than 15*p

		for (size_t i = static_sieve_primes.size() + 1; small_primes_lookup[i] <= last_prime_for_stepping_by_fifteen; ++i)
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

		for (int i = 0; i < 8; ++i)
		{
			masks[i].set_all();

			masks[i].clear_bit(i + (p - p));
			masks[i].clear_bit(i + (2 * p) - p);
			masks[i].clear_bit(i + (4 * p) - p);
			masks[i].clear_bit(i + (7 * p) - p);
			masks[i].clear_bit(i + (8 * p) - p);

			if constexpr (p <= 23)
			{
				masks[i].clear_bit(i + (11 * p) - p);
			}

			if constexpr (p <= 19)
			{
				masks[i].clear_bit(i + (13 * p) - p);
				masks[i].clear_bit(i + (14 * p) - p);
			}
		}

		return masks;
	}
	static constexpr std::array<bit_array<256>, 8> sieve_masks_p19 = generate_sieve_masks<19>();
	static constexpr std::array<bit_array<256>, 8> sieve_masks_p23 = generate_sieve_masks<23>();
	static constexpr std::array<bit_array<256>, 8> sieve_masks_p29 = generate_sieve_masks<29>();
	static constexpr std::array<bit_array<256>, 8> sieve_masks_p31 = generate_sieve_masks<31>();

	template<size_t p>
	consteval std::array<bit_array<256>, 16> generate_wide_sieve_masks()
	{
		decltype(generate_wide_sieve_masks<p>()) masks{};

		const size_t mask_1_offset = 1 * p;
		const size_t mask_2_offset = 8 * p;

		for (size_t j = 0; j < 2ull * 8; j += 2)
		{
			// Create a 256-bit mask to mask against 1p, 2p, 4p, and possibly 7p.

			masks[j].set_all();
			masks[j].clear_bit((j / 2) + (1 * p) - mask_1_offset);
			masks[j].clear_bit((j / 2) + (2 * p) - mask_1_offset);
			masks[j].clear_bit((j / 2) + (4 * p) - mask_1_offset);
			if constexpr (p == 37 || p == 41)
				masks[j].clear_bit((j / 2) + (7 * p) - mask_1_offset);

			// Create a second mask to mask against 8p, 11p, 13p, and possibly 14p.

			// The bit offset for 1p is j/2. Calculate the bit offset for 8p
			const size_t bit_off = ((j / 2) + mask_2_offset - mask_1_offset) % 8;

			masks[j + 1].set_all();
			masks[j + 1].clear_bit(bit_off + (8 * p) - mask_2_offset);
			masks[j + 1].clear_bit(bit_off + (11 * p) - mask_2_offset);
			masks[j + 1].clear_bit(bit_off + (13 * p) - mask_2_offset);
			if constexpr (p == 37 || p == 41)
				masks[j + 1].clear_bit(bit_off + (14 * p) - mask_2_offset);
		}

		return masks;
	}
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p37 = generate_wide_sieve_masks<37>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p41 = generate_wide_sieve_masks<41>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p43 = generate_wide_sieve_masks<43>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p47 = generate_wide_sieve_masks<47>();

	template<size_t p>
	consteval std::array<bit_array<256>, 16> generate_wider_sieve_masks()
	{
		decltype(generate_wider_sieve_masks<p>()) masks{};

		const size_t mask_1_offset = 1 * p;
		const size_t mask_2_offset = 7 * p;

		for (size_t j = 0; j < 2ull * 8; j += 2)
		{
			// Create a 256-bit mask to mask against 1p, 2p, and 4p.

			masks[j].set_all();
			masks[j].clear_bit((j / 2) + (1 * p) - mask_1_offset);
			masks[j].clear_bit((j / 2) + (2 * p) - mask_1_offset);
			masks[j].clear_bit((j / 2) + (4 * p) - mask_1_offset);

			// Create a second mask to mask against 7p, 8p, and 11p.

			// The bit offset for 1p is j/2. Calculate the bit offset for 7p
			const size_t bit_off = ((j / 2) + mask_2_offset - mask_1_offset) % 8;

			masks[j + 1].set_all();
			masks[j + 1].clear_bit(bit_off + (7 * p) - mask_2_offset);
			masks[j + 1].clear_bit(bit_off + (8 * p) - mask_2_offset);
			masks[j + 1].clear_bit(bit_off + (11 * p) - mask_2_offset);
		}

		return masks;
	}
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p53 = generate_wider_sieve_masks<53>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p59 = generate_wider_sieve_masks<59>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p61 = generate_wider_sieve_masks<61>();

	template<size_t p>
	consteval std::array<bit_array<256>, 16> generate_even_wider_sieve_masks()
	{
		decltype(generate_wider_sieve_masks<p>()) masks{};

		const size_t mask_1_offset = 1 * p;
		const size_t mask_2_offset = 11 * p;

		for (size_t j = 0; j < 2ull * 8; j += 2)
		{
			// Create a 256-bit mask to mask against 1p, 2p, and 4p.

			masks[j].set_all();
			masks[j].clear_bit((j / 2) + (1 * p) - mask_1_offset);
			masks[j].clear_bit((j / 2) + (2 * p) - mask_1_offset);
			masks[j].clear_bit((j / 2) + (4 * p) - mask_1_offset);

			// Create a second mask to mask against 11p, 13p, and 14p.

			// The bit offset for 1p is j/2. Calculate the bit offset for 7p
			const size_t bit_off = ((j / 2) + mask_2_offset - mask_1_offset) % 8;

			masks[j + 1].set_all();
			masks[j + 1].clear_bit(bit_off + (11 * p) - mask_2_offset);
			masks[j + 1].clear_bit(bit_off + (13 * p) - mask_2_offset);
			masks[j + 1].clear_bit(bit_off + (14 * p) - mask_2_offset);
		}

		return masks;
	}
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p67 = generate_even_wider_sieve_masks<67>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p71 = generate_even_wider_sieve_masks<71>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p73 = generate_even_wider_sieve_masks<73>();
	static constexpr std::array<bit_array<256>, 16> sieve_masks_p79 = generate_even_wider_sieve_masks<79>();

	template<size_t p>
	__forceinline void vectorized_sieve_pass(sieve_container& sieve,
											 const std::array<bit_array<256>, 8>& sieve_masks,
											 const sieve_prime_t*& prime_ptr,
											 sieve_offset_t*& offset_cache_ptr)
	{
		constexpr size_t sieve_end = sieve_container::size();
		constexpr size_t padded_end = sieve_end - (15 * p);

		// Get the position of the next odd multiple of p*15
		size_t j = *offset_cache_ptr;

		do
		{
			const size_t byte_index = (j + p) / 8;
			const size_t bit_index = (j + p) % 8; // 0..7

			uint256_t* const sieve_ptr = (uint256_t*)(sieve.data() + byte_index);

			const uint256_t mask = _mm256_loadu_si256((const uint256_t*)sieve_masks[bit_index].data());

			uint256_t sieve_data = _mm256_loadu_si256(sieve_ptr);
			sieve_data = _mm256_and_si256(mask, sieve_data);
			_mm256_storeu_si256(sieve_ptr, sieve_data);

			if constexpr (p >= 29)
			{
				sieve.clear_bit(j + 11 * p);
			}

			if constexpr (p >= 23)
			{
				sieve.clear_bit(j + 13 * p);
				sieve.clear_bit(j + 14 * p);
			}

			j += (15 * p);
			// Stop marking 15*p early (don't handle padding)
		} while (j < padded_end);

		// Calculate and cache the offset for the next sieving
		*offset_cache_ptr = sieve_offset_t((j + (15 * p)) - sieve_end);

		++prime_ptr;
		++offset_cache_ptr;
	}

	template<size_t p>
	__forceinline void wide_vectorized_sieve_pass(sieve_container& sieve,
												  const std::array<bit_array<256>, 16>& sieve_masks,
												  const sieve_prime_t*& prime_ptr,
												  sieve_offset_t*& offset_cache_ptr)
	{
		constexpr size_t sieve_end = sieve_container::size();
		constexpr size_t padded_end = sieve_end - (15 * p);

		// Get the position of the next odd multiple of p*15
		size_t j = *offset_cache_ptr;

		do
		{
			const size_t bit_index = (j + p) % 8; // 0..7
			const uint256_t mask_lo = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index * 2].data());
			const uint256_t mask_hi = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index * 2 + 1].data());

			const size_t byte_index_lo = (j + p) / 8;
			const size_t byte_index_hi = (j + (8 * p)) / 8;
			uint256_t* const sieve_ptr_lo = (uint256_t*)(sieve.data() + byte_index_lo);
			uint256_t* const sieve_ptr_hi = (uint256_t*)(sieve.data() + byte_index_hi);

			uint256_t sieve_data_lo = _mm256_loadu_si256(sieve_ptr_lo);
			sieve_data_lo = _mm256_and_si256(mask_lo, sieve_data_lo);
			_mm256_storeu_si256(sieve_ptr_lo, sieve_data_lo);

			uint256_t sieve_data_hi = _mm256_loadu_si256(sieve_ptr_hi);
			sieve_data_hi = _mm256_and_si256(mask_hi, sieve_data_hi);
			_mm256_storeu_si256(sieve_ptr_hi, sieve_data_hi);

			if constexpr (p == 43 || p == 47)
			{
				sieve.clear_bit(j + (7 * p));
				sieve.clear_bit(j + (14 * p));
			}

			j += (15 * p);
		} while (j < padded_end);

		// Calculate and cache the offset for the next sieving
		*offset_cache_ptr = sieve_offset_t((j + (15 * p)) - sieve_end);

		++prime_ptr;
		++offset_cache_ptr;
	}

	template<size_t p>
	__forceinline void wider_vectorized_sieve_pass(sieve_container& sieve,
												   const std::array<bit_array<256>, 16>& sieve_masks,
												   const sieve_prime_t*& prime_ptr,
												   sieve_offset_t*& offset_cache_ptr)
	{
		constexpr size_t sieve_end = sieve_container::size();
		constexpr size_t padded_end = sieve_end - (15 * p);

		// Get the position of the next odd multiple of p*15
		size_t j = *offset_cache_ptr;

		do
		{
			const size_t bit_index = (j + p) % 8; // 0..7
			const uint256_t mask_lo = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index * 2].data());
			const uint256_t mask_hi = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index * 2 + 1].data());

			const size_t byte_index_lo = (j + p) / 8;
			const size_t byte_index_hi = (j + (7 * p)) / 8;
			uint256_t* const sieve_ptr_lo = (uint256_t*)(sieve.data() + byte_index_lo);
			uint256_t* const sieve_ptr_hi = (uint256_t*)(sieve.data() + byte_index_hi);

			// mask against 1*p, 2*p, and 4*p
			uint256_t sieve_data_lo = _mm256_loadu_si256(sieve_ptr_lo);
			sieve_data_lo = _mm256_and_si256(mask_lo, sieve_data_lo);
			_mm256_storeu_si256(sieve_ptr_lo, sieve_data_lo);

			// mask against 7*p, 8*p, 11*p
			uint256_t sieve_data_hi = _mm256_loadu_si256(sieve_ptr_hi);
			sieve_data_hi = _mm256_and_si256(mask_hi, sieve_data_hi);
			_mm256_storeu_si256(sieve_ptr_hi, sieve_data_hi);

			sieve.clear_bit(j + (13 * p));
			sieve.clear_bit(j + (14 * p));

			j += (15 * p);
		} while (j < padded_end);

		// Calculate and cache the offset for the next sieving
		*offset_cache_ptr = sieve_offset_t((j + (15 * p)) - sieve_end);

		++prime_ptr;
		++offset_cache_ptr;
	}

	template<size_t p>
	__forceinline void even_wider_vectorized_sieve_pass(sieve_container& sieve,
														const std::array<bit_array<256>, 16>& sieve_masks,
														const sieve_prime_t*& prime_ptr,
														sieve_offset_t*& offset_cache_ptr)
	{
		constexpr size_t sieve_end = sieve_container::size();
		constexpr size_t padded_end = sieve_end - (15 * p);

		// Get the position of the next odd multiple of p*15
		size_t j = *offset_cache_ptr;

		do
		{
			const size_t bit_index = (j + p) % 8; // 0..7
			const uint256_t mask_lo = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index * 2].data());
			const uint256_t mask_hi = _mm256_loadu_si256((uint256_t*)sieve_masks[bit_index * 2 + 1].data());

			const size_t byte_index_lo = (j + p) / 8;
			const size_t byte_index_hi = (j + (11 * p)) / 8;
			uint256_t* const sieve_ptr_lo = (uint256_t*)(sieve.data() + byte_index_lo);
			uint256_t* const sieve_ptr_hi = (uint256_t*)(sieve.data() + byte_index_hi);

			// mask against 1*p, 2*p, and 4*p
			uint256_t sieve_data_lo = _mm256_loadu_si256(sieve_ptr_lo);
			sieve_data_lo = _mm256_and_si256(mask_lo, sieve_data_lo);
			_mm256_storeu_si256(sieve_ptr_lo, sieve_data_lo);

			// mask against 11*p, 13*p, 14*p
			uint256_t sieve_data_hi = _mm256_loadu_si256(sieve_ptr_hi);
			sieve_data_hi = _mm256_and_si256(mask_hi, sieve_data_hi);
			_mm256_storeu_si256(sieve_ptr_hi, sieve_data_hi);

			sieve.clear_bit(j + (7 * p));
			sieve.clear_bit(j + (8 * p));

			j += (15 * p);
		} while (j < padded_end);

		// Calculate and cache the offset for the next sieving
		*offset_cache_ptr = sieve_offset_t((j + (15 * p)) - sieve_end);

		++prime_ptr;
		++offset_cache_ptr;
	}

	void partial_sieve(const uint64_t number,
					   sieve_container& sieve
					   count_passes(, size_t& ps15))
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

		if constexpr (small_primes_lookup[static_sieve_primes.size() + 1] == 19)
		{
			// one simd write, plus 0-3 scalar writes for 11,13,14
			vectorized_sieve_pass<19>(sieve, sieve_masks_p19, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<23>(sieve, sieve_masks_p23, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<29>(sieve, sieve_masks_p29, prime_ptr, offset_cache_ptr);
			vectorized_sieve_pass<31>(sieve, sieve_masks_p31, prime_ptr, offset_cache_ptr);

			// two simd writes, plus 0-2 scalar writes for 7,14
			wide_vectorized_sieve_pass<37>(sieve, sieve_masks_p37, prime_ptr, offset_cache_ptr);
			wide_vectorized_sieve_pass<41>(sieve, sieve_masks_p41, prime_ptr, offset_cache_ptr);
			wide_vectorized_sieve_pass<43>(sieve, sieve_masks_p43, prime_ptr, offset_cache_ptr);
			wide_vectorized_sieve_pass<47>(sieve, sieve_masks_p47, prime_ptr, offset_cache_ptr);

			// two simd writes, plus 2 scalar writes for 13,14
			wider_vectorized_sieve_pass<53>(sieve, sieve_masks_p53, prime_ptr, offset_cache_ptr);
			wider_vectorized_sieve_pass<59>(sieve, sieve_masks_p59, prime_ptr, offset_cache_ptr);
			wider_vectorized_sieve_pass<61>(sieve, sieve_masks_p61, prime_ptr, offset_cache_ptr);

			// two simd writes, plus 2 scalar writes for 7,8
			even_wider_vectorized_sieve_pass<67>(sieve, sieve_masks_p67, prime_ptr, offset_cache_ptr);
			even_wider_vectorized_sieve_pass<71>(sieve, sieve_masks_p71, prime_ptr, offset_cache_ptr);
			even_wider_vectorized_sieve_pass<73>(sieve, sieve_masks_p73, prime_ptr, offset_cache_ptr);
			even_wider_vectorized_sieve_pass<79>(sieve, sieve_masks_p79, prime_ptr, offset_cache_ptr);
		}

		double density = double(sieve.count_bits()) / sieve_container::size();

		for (;;)
		{
			// If we stop sieving early, we still need to update our offsets cache
			if (density < density_threshold)
			{
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

			// Each iteration, j points to an (already marked) multiple of p*15.
			// Mark off multiples of p.
			do
			{
				sieve.clear_bit(j + p);
				sieve.clear_bit(j + 2 * p);
				sieve.clear_bit(j + 4 * p);
				sieve.clear_bit(j + 7 * p);
				sieve.clear_bit(j + 8 * p);
				sieve.clear_bit(j + 11 * p);
				sieve.clear_bit(j + 13 * p);
				sieve.clear_bit(j + 14 * p);

				j += (15 * p);
			} while (j < padded_end);

			// Calculate and cache the offset for the next sieving
			*offset_cache_ptr++ = sieve_offset_t((j + (15 * p)) - sieve_end);

			if (p == last_prime_for_stepping_by_fifteen) break;
		}

		count_passes(ps15 += sieve.count_bits());
	}



	namespace detail
	{
		__forceinline void gather_bits_step(size_t& bit_idx, uint64_t& chunk,
											size_t*& candidates, const size_t number,
											const sieve_container& sieve)
		{
			const size_t low_bits = bit_idx & 0b111;

			// Zero out the bottom n bits, where n is our offset into the byte
			const uint64_t masked_chunk = chunk & (uint64_t(-1) << low_bits);
			const size_t trailing_zeroes = _tzcnt_u64(masked_chunk);

			bit_idx -= low_bits;
			bit_idx += trailing_zeroes;

			*candidates = number + (2 * bit_idx); // always store to avoid a branch

			const size_t chunk_is_non_zero = (trailing_zeroes >> 6) ^ 1; // max tzc == 64 == 0b0100'0000. Shift by 6 then invert.

			candidates += chunk_is_non_zero; // if we found a bit, keep the candidate 
			bit_idx += chunk_is_non_zero; // if we found a bit, our next index is one past it

			chunk = sieve.load_u64_chunk(bit_idx); // load ahead
		}

		__forceinline void gather_bits_cleanup(size_t& bit_idx, const size_t block_end_idx, uint64_t& chunk,
											   size_t*& candidates, const size_t number)
		{
			const size_t low_bits = bit_idx & 0b111;
			chunk >>= low_bits; // offset range 0-7

			// we have 0-63 bits left - mask out any other (high) bits
			chunk &= (1ull << (block_end_idx - bit_idx)) - 1;

			while (chunk != 0)
			{
				const size_t trailing_zeroes = _tzcnt_u64(chunk);
				chunk >>= trailing_zeroes; // always 0-63
				chunk >>= 1; // one more shift to clear the bit

				bit_idx += trailing_zeroes;

				*candidates = number + (2 * bit_idx); // always store to avoid a branch
				++candidates; // always increment (since chunk != 0, we have a candidate)

				++bit_idx; // step past current candidate
			}
		}
	}

	size_t* gather_sieve_results(size_t* candidates,
								 const sieve_container& sieve,
								 const size_t number)
	{
		// one block contains many 64-bit chunks
		constexpr size_t block_size = ((sieve_container::size() / 4) / 8) * 8;
		constexpr size_t chunk_size = 64;

		size_t block_0_end = block_size * 1;
		size_t block_1_end = block_size * 2;
		size_t block_2_end = block_size * 3;
		size_t block_3_end = sieve.size();

		size_t block_0_idx = 0;
		size_t block_1_idx = block_size * 1;
		size_t block_2_idx = block_size * 2;
		size_t block_3_idx = block_size * 3;

		uint64_t chunk_0 = sieve.load_u64_chunk(block_0_idx);
		uint64_t chunk_1 = sieve.load_u64_chunk(block_1_idx);
		uint64_t chunk_2 = sieve.load_u64_chunk(block_2_idx);
		uint64_t chunk_3 = sieve.load_u64_chunk(block_3_idx);

		// Gather from 4 locations in parallel

		size_t safe_reads = block_size / chunk_size;

		do
		{
			do
			{
				detail::gather_bits_step(block_0_idx, chunk_0, candidates, number, sieve);
				detail::gather_bits_step(block_1_idx, chunk_1, candidates, number, sieve);
				detail::gather_bits_step(block_2_idx, chunk_2, candidates, number, sieve);
				detail::gather_bits_step(block_3_idx, chunk_3, candidates, number, sieve);
			} while (--safe_reads != 0);

			safe_reads = util::min(block_0_end - block_0_idx,
								   block_1_end - block_1_idx,
								   block_2_end - block_2_idx,
								   block_3_end - block_3_idx) / chunk_size;
		} while (safe_reads != 0);

		// At least one idx is within 64 bits of the end of its block.
		// If it is one of the first three blocks, swap with block 3 for cleanup

		if (block_0_end - block_0_idx < chunk_size)
		{
			std::swap(block_0_idx, block_3_idx);
			std::swap(block_0_end, block_3_end);
			std::swap(chunk_0, chunk_3);
		}
		else if (block_1_end - block_1_idx < chunk_size)
		{
			std::swap(block_1_idx, block_3_idx);
			std::swap(block_1_end, block_3_end);
			std::swap(chunk_1, chunk_3);
		}
		else if (block_2_end - block_2_idx < chunk_size)
		{
			std::swap(block_2_idx, block_3_idx);
			std::swap(block_2_end, block_3_end);
			std::swap(chunk_2, chunk_3);
		}

		// in all cases, finish whatever block 3 is now
		detail::gather_bits_cleanup(block_3_idx, block_3_end, chunk_3, candidates, number);



		// Gather from 3 locations in parallel

		// recalculate
		safe_reads = util::min(block_0_end - block_0_idx,
							   block_1_end - block_1_idx,
							   block_2_end - block_2_idx) / chunk_size;

		while (safe_reads != 0)
		{
			do
			{
				detail::gather_bits_step(block_0_idx, chunk_0, candidates, number, sieve);
				detail::gather_bits_step(block_1_idx, chunk_1, candidates, number, sieve);
				detail::gather_bits_step(block_2_idx, chunk_2, candidates, number, sieve);
			} while (--safe_reads != 0);

			safe_reads = util::min(block_0_end - block_0_idx,
								   block_1_end - block_1_idx,
								   block_2_end - block_2_idx) / chunk_size;
		}

		if (block_0_end - block_0_idx < chunk_size)
		{
			std::swap(block_0_idx, block_2_idx);
			std::swap(block_0_end, block_2_end);
			std::swap(chunk_0, chunk_2);
		}
		else if (block_1_end - block_1_idx < chunk_size)
		{
			std::swap(block_1_idx, block_2_idx);
			std::swap(block_1_end, block_2_end);
			std::swap(chunk_1, chunk_2);
		}

		detail::gather_bits_cleanup(block_2_idx, block_2_end, chunk_2, candidates, number);



		// Gather from 2 locations in parallel

		safe_reads = util::min(block_0_end - block_0_idx,
							   block_1_end - block_1_idx) / chunk_size;

		while (safe_reads != 0)
		{
			do
			{
				detail::gather_bits_step(block_0_idx, chunk_0, candidates, number, sieve);
				detail::gather_bits_step(block_1_idx, chunk_1, candidates, number, sieve);
			} while (--safe_reads != 0);

			safe_reads = util::min(block_0_end - block_0_idx,
								   block_1_end - block_1_idx) / chunk_size;
		}

		if (block_0_end - block_0_idx < chunk_size)
		{
			std::swap(block_0_idx, block_1_idx);
			std::swap(block_0_end, block_1_end);
			std::swap(chunk_0, chunk_1);
		}

		detail::gather_bits_cleanup(block_1_idx, block_1_end, chunk_1, candidates, number);



		// Gather from remaining idx

		safe_reads = (block_0_end - block_0_idx) / chunk_size;

		while (safe_reads != 0)
		{
			do
			{
				detail::gather_bits_step(block_0_idx, chunk_0, candidates, number, sieve);
			} while (--safe_reads != 0);

			safe_reads = (block_0_end - block_0_idx) / chunk_size;
		}

		detail::gather_bits_cleanup(block_0_idx, block_0_end, chunk_0, candidates, number);



		return candidates;
	}

}
