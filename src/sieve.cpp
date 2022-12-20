
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

	void verify_sieve_offset_cache(const uint64_t start)
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

			if (p == largest_sieve_prime) break;
		}

		count_passes(ps15 += sieve.count_bits());
	}



	namespace detail
	{
		constexpr uint64_t n_sieve_chunks = (sieve_container::size() / 64) + (sieve_container::size() % 64 != 0);
		constexpr size_t max_pc = 27; // 27 is the highest popcount of 64-bit chunks in the static sieve

		using chunk_count_t = util::narrowest_uint_for_val<n_sieve_chunks>;
		static std::array<chunk_count_t, max_pc + 1> n_chunks_with_pc{};
		// [popcount][chunk]
		static std::array<std::array<uint64_t, n_sieve_chunks>, max_pc + 1> sorted_chunks{};

		using chunk_idx_t = util::narrowest_uint_for_val<n_sieve_chunks>;
		// [popcount][idx]
		static std::array<std::array<chunk_idx_t, n_sieve_chunks>, max_pc + 1> chunk_indexes{};

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
