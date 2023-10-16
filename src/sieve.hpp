#pragma once

#include <iostream>
#include <iomanip>

#include "math/math.hpp"
#include "util/types.hpp"

namespace mbp::prime_sieve
{
	const sieve_container generate_static_sieve();

	using sieve_offset_t = util::narrowest_uint_for_val<static_sieve_size>;
	extern std::array<sieve_offset_t, small_primes_lookup.size()> sieve_offsets_cache;

	void set_up_sieve_offsets_cache(const size_t start);

	// precalculate stride - (sieve_size % stride)
	constexpr static std::array<sieve_offset_t, small_primes_lookup.size()> modulo_precomp = []() consteval {
		std::array<sieve_offset_t, small_primes_lookup.size()> arr{};

		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			const size_t prime = small_primes_lookup[i];

			if (prime >= 17 && prime <= largest_aligned_vector_sieve_prime) // handle separately
			{
				arr[i] = sieve_offset_t(prime - (sieve_container::size() % prime));
				continue;
			}

			size_t stride = prime * 15;
			if (prime <= largest_vector_sieve_prime) stride *= 8;
			arr[i] = sieve_offset_t(stride - (sieve_container::size() % stride));
		}

		return arr;
	}();

	static void update_sieve_offsets_cache(const sieve_prime_t* prime_ptr,
										   sieve_offset_t* offset_ptr)
	{
		const auto* mp_ptr = modulo_precomp.data() + (prime_ptr - small_primes_lookup.data());

		//std::cout << __FILE__ << ' ' << __LINE__ << '\n';
		//std::cout << "small_primes_lookup.data() = " << small_primes_lookup.data() << '\n';
		//std::cout << "prime_ptr =                  " << prime_ptr << '\n';
		//std::cout << "difference = " << prime_ptr - small_primes_lookup.data() << '\n';
		//std::cin.ignore();

		// Calculate the offset of the next odd multiple of the stride size.
		for (size_t prime = *prime_ptr;
			 prime <= largest_sieve_prime;
			 prime = *++prime_ptr, ++offset_ptr, ++mp_ptr)
		{
			size_t stride = prime;

			// *15 if the prime is used for unaligned vector sieving or scalar sieving
			stride *= (prime > largest_aligned_vector_sieve_prime) ? 15 : 1;

			// *8 if the prime is used for unaligned vector sieving
			stride *= (prime > largest_aligned_vector_sieve_prime &&
					   prime <= largest_vector_sieve_prime) ? 8 : 1;

			size_t n = *offset_ptr;
			n += *mp_ptr;
			n = util::min(n, n - stride);

			*offset_ptr = sieve_offset_t(n);
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
	consteval std::array<bit_array<256>, p> generate_aligned_sieve_masks()
	{
		std::array<bit_array<256>, p> masks{};

		for (size_t bit_offset = 0; bit_offset < p; ++bit_offset)
		{
			masks[bit_offset].set_all();

			for (size_t i = bit_offset; i < 256; i += p)
			{
				masks[bit_offset].clear_bit(i);
			}
		}

		return masks;
	}



	template<size_t p, size_t chunk = 0>
	__forceinline void generate_aligned_vector_writes(
		uint256_t* const ptr, const auto& masks,
		const uint256_t& mask_0, const uint256_t& mask_1, const uint256_t& mask_2, const uint256_t& mask_3,
		const uint256_t& mask_4, const uint256_t& mask_5, const uint256_t& mask_6, const uint256_t& mask_7,
		const uint256_t& mask_8, const uint256_t& mask_9, const uint256_t& mask_10, const uint256_t& mask_11,
		const uint256_t& mask_12, const uint256_t& mask_13)
	{
		constexpr size_t bit_offset = (p - ((chunk * 256) % p)) % p; // cleaner way to do this?

		uint256_t mask{};

		if constexpr (bit_offset == 0) mask = mask_0;
		else if constexpr (bit_offset == 1) mask = mask_1;
		else if constexpr (bit_offset == 2) mask = mask_2;
		else if constexpr (bit_offset == 3) mask = mask_3;
		else if constexpr (bit_offset == 4) mask = mask_4;
		else if constexpr (bit_offset == 5) mask = mask_5;
		else if constexpr (bit_offset == 6) mask = mask_6;
		else if constexpr (bit_offset == 7) mask = mask_7;
		else if constexpr (bit_offset == 8) mask = mask_8;
		else if constexpr (bit_offset == 9) mask = mask_9;
		else if constexpr (bit_offset == 10) mask = mask_10;
		else if constexpr (bit_offset == 11) mask = mask_11;
		else if constexpr (bit_offset == 12) mask = mask_12;
		else if constexpr (bit_offset == 13) mask = mask_13;
		else mask = *(uint256_t*)masks[bit_offset].data();

		*ptr = _mm256_and_si256(mask, *ptr);

		if constexpr (chunk + 1 < p)
			generate_aligned_vector_writes<p, chunk + 1>(
				ptr + 1, masks,
				mask_0, mask_1, mask_2, mask_3, mask_4,
				mask_5, mask_6, mask_7, mask_8, mask_9,
				mask_10, mask_11, mask_12, mask_13);
	}

	template<size_t p>
	__forceinline void aligned_vectorized_sieve_pass(sieve_container& sieve,
													 const sieve_prime_t*& prime_ptr,
													 sieve_offset_t*& offset_cache_ptr)
	{
		constexpr static std::array masks = generate_aligned_sieve_masks<p>();

		size_t j = *offset_cache_ptr;
		uint256_t* ptr = (uint256_t*)sieve.data();

		// align j with a bit offset of 0
		while (j != 0)
		{
			*ptr = _mm256_and_si256(*(uint256_t*)masks[j].data(), *ptr);

			j += p - (256 % p);
			j = util::min(j, j - p);
			++ptr;
		}

		// MSVC will use 14 of 16 YMM registers for constants without stack spilling
		const uint256_t mask_0 = _mm256_loadu_si256((uint256_t*)masks[0].data());
		const uint256_t mask_1 = _mm256_loadu_si256((uint256_t*)masks[1].data());
		const uint256_t mask_2 = _mm256_loadu_si256((uint256_t*)masks[2].data());
		const uint256_t mask_3 = _mm256_loadu_si256((uint256_t*)masks[3].data());
		const uint256_t mask_4 = _mm256_loadu_si256((uint256_t*)masks[4].data());
		const uint256_t mask_5 = _mm256_loadu_si256((uint256_t*)masks[5].data());
		const uint256_t mask_6 = _mm256_loadu_si256((uint256_t*)masks[6].data());
		const uint256_t mask_7 = _mm256_loadu_si256((uint256_t*)masks[7].data());
		const uint256_t mask_8 = _mm256_loadu_si256((uint256_t*)masks[8].data());
		const uint256_t mask_9 = _mm256_loadu_si256((uint256_t*)masks[9].data());
		const uint256_t mask_10 = _mm256_loadu_si256((uint256_t*)masks[10].data());
		const uint256_t mask_11 = _mm256_loadu_si256((uint256_t*)masks[11].data());
		const uint256_t mask_12 = _mm256_loadu_si256((uint256_t*)masks[12].data());
		const uint256_t mask_13 = _mm256_loadu_si256((uint256_t*)masks[13].data());

		// iterate until we reach the last p-1 chunks
		constexpr size_t n_chunks = sieve_container::size() / 256ull;
		const uint256_t* const aligned_end = ((uint256_t*)sieve.data()) + n_chunks;
		const uint256_t* const extra_aligned_end = aligned_end - (p - 1);

		do
		{
			generate_aligned_vector_writes<p>(ptr, masks,
											  mask_0, mask_1, mask_2, mask_3, mask_4,
											  mask_5, mask_6, mask_7, mask_8, mask_9,
											  mask_10, mask_11, mask_12, mask_13);
			ptr += p; // advance by p*32 bytes
		} while (ptr < extra_aligned_end);

		while (ptr < aligned_end)
		{
			*ptr = _mm256_and_si256(*(uint256_t*)masks[j].data(), *ptr);

			j += p - (256 % p);
			j = util::min(j, j - p);
			++ptr;
		}

		j += p - ((sieve_container::size() % 256) % p);
		j = util::min(j, j - p);

		*offset_cache_ptr = sieve_offset_t(j);

		++prime_ptr;
		++offset_cache_ptr;
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
								   const uint8_t mask)
	{
		size_t next_offset = offset + delta_to_next_offset;

		offset >>= 3;
		sieve_data[offset] &= mask;

		return next_offset;
	}

	__forceinline size_t clear_8_of_15(const size_t p, size_t j, uint8_t* const sieve_data,
									   const uint8_t mask_a,
									   const uint8_t mask_b,
									   const uint8_t mask_c,
									   const uint8_t mask_d,
									   const uint8_t mask_e,
									   const uint8_t mask_f,
									   const uint8_t mask_g,
									   const uint8_t mask_h)
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

	template<size_t stride, size_t m_offset>
	__forceinline void generate_1_scalar_write(const size_t p, size_t& byte_offset, uint8_t* const p1_ptr)
	{
		constexpr size_t p_mod_8 = 1; // fixme

		// Calculate the adjustment needed to get the basic byte_offset to the actual byte
		constexpr size_t adjust = ((stride * 15) + m_offset - 1) / 8;

		constexpr size_t element_offset = (p_mod_8 * 15 * stride) + (p_mod_8 * (m_offset - 1));
		constexpr uint8_t mask = ~(1u << (element_offset % 8)) & 0xFF;

		p1_ptr[byte_offset + adjust] &= mask;

		const size_t step = p / 8; // calculated once up front

		// replace this by calculating the difference between the current and next m_offsets
		if constexpr (m_offset == 1) byte_offset += step;
		if constexpr (m_offset == 2) byte_offset += step * 2;
		if constexpr (m_offset == 4) byte_offset += step * 3;
		if constexpr (m_offset == 7) byte_offset += step;
		if constexpr (m_offset == 8) byte_offset += step * 3;
		if constexpr (m_offset == 11) byte_offset += step * 2;
		if constexpr (m_offset == 13) byte_offset += step;
		if constexpr (m_offset == 14) byte_offset += step * 2;
	}

	// for now, assume p % 8 == 1
	__forceinline void generate_64_scalar_writes(const size_t p, uint8_t* const p1_ptr)
	{
		constexpr size_t m_offsets[]{ 1, 2, 4, 7, 8, 11, 13, 14 };

		size_t b = 0;

		generate_1_scalar_write<0, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<0, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<1, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<1, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<2, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<2, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<3, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<3, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<4, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<4, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<5, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<5, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<6, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<6, m_offsets[7]>(p, b, p1_ptr);

		generate_1_scalar_write<7, m_offsets[0]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[1]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[2]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[3]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[4]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[5]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[6]>(p, b, p1_ptr);
		generate_1_scalar_write<7, m_offsets[7]>(p, b, p1_ptr);
	}



	inline_toggle static void partial_sieve(sieve_container& sieve,
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

		// don't do any sieving if our bitmasks + static sieve already cleared enough
		if (density < vector_density_threshold)
		{
			//std::cout << __FILE__ << ' ' << __LINE__ << '\n';
			//std::cout << "small_primes_lookup.data() = " << small_primes_lookup.data() << '\n';
			//std::cout << "prime_ptr =                  " << prime_ptr << '\n';
			//std::cout << "difference = " << prime_ptr - small_primes_lookup.data() << '\n';
			//std::cin.ignore();

			update_sieve_offsets_cache(prime_ptr, offset_cache_ptr);
			return;
		}

		//if (*offset_cache_ptr > 20'000)
		//{
		//	std::cout << "Hol' up. Dumping sieve_offsets_cache:\n";

		//	for (size_t i = 0; i < sieve_offsets_cache.size(); ++i)
		//	{
		//		std::cout << i << '\t' << sieve_offsets_cache[i] << '\n';
		//	}

		//	std::cin.ignore();
		//}

		static_assert(small_primes_lookup[static_sieve_primes.size() + 1] == 17);
		aligned_vectorized_sieve_pass<17>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<19>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<23>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<29>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<31>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<37>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<41>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<43>(sieve, prime_ptr, offset_cache_ptr);
		aligned_vectorized_sieve_pass<47>(sieve, prime_ptr, offset_cache_ptr);
		static_assert(largest_aligned_vector_sieve_prime == 47);
		vectorized_sieve_pass<53>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<59>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<61>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<67>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<71>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<73>(sieve, prime_ptr, offset_cache_ptr);
		vectorized_sieve_pass<79>(sieve, prime_ptr, offset_cache_ptr);
		static_assert(largest_vector_sieve_prime == 79);

		constexpr double vectorized_sieving_removes = .3106; // 31.06% (may be out of date)
		constexpr double scale = 1.0 - vectorized_sieving_removes;

		density *= scale;

		for (;;)
		{
			if (density < scalar_density_threshold)
			{
				//std::cout << __FILE__ << ' ' << __LINE__ << '\n';
				//std::cout << "small_primes_lookup.data() = " << small_primes_lookup.data() << '\n';
				//std::cout << "prime_ptr =                  " << prime_ptr << '\n';
				//std::cout << "difference = " << prime_ptr - small_primes_lookup.data() << '\n';
				//std::cin.ignore();

				// If we stop sieving early, we still need to update our offsets cache
				update_sieve_offsets_cache(prime_ptr, offset_cache_ptr);
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

			if (false) // if (p % 8 == 1)
			{
				do
				{
					generate_64_scalar_writes(p, sieve_data + ((j + p) / 8));

					j += 15ull * 8 * p;
				} while (j < extra_padded_end);
			}
			else // p % 8 == 3, 5, or 7
			{
				// generate eight byte masks in usage order, starting from (j + p) % 8 == 0
				const uint8_t mask_a = ~(1ull << ((j + (1ull * p)) % 8));
				const uint8_t mask_b = ~(1ull << ((j + (2ull * p)) % 8));
				const uint8_t mask_c = ~(1ull << ((j + (4ull * p)) % 8));
				const uint8_t mask_d = ~(1ull << ((j + (7ull * p)) % 8));
				const uint8_t mask_e = ~(1ull << ((j + (8ull * p)) % 8));
				const uint8_t mask_f = ~(1ull << ((j + (11ull * p)) % 8));
				const uint8_t mask_g = ~(1ull << ((j + (13ull * p)) % 8));
				const uint8_t mask_h = ~(1ull << ((j + (14ull * p)) % 8));

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
			}

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

		using chunk_count_t = util::narrowest_uint_for_val<n_sieve_chunks>;
		using chunk_idx_t = chunk_count_t;

		// [popcount] -> number of chunks with that popcount
		static std::array<chunk_count_t, 8> n_chunks_with_pc;
		// [popcount][chunks]
		static std::array<std::array<uint64_t, n_sieve_chunks>, 8> sorted_chunks;
		// [popcount][chunk indexes]
		static std::array<std::array<chunk_idx_t, n_sieve_chunks>, 8> chunk_indexes;

		static const std::array<uint64_t, 16> permute_masks = [] {
			std::array<uint64_t, 16> masks{};

			for (size_t i = 0; i < 16; ++i)
			{
				const uint64_t byte_mask = _pdep_u64(~i, 0x0001000100010001) * 0xFFFF; // convert bits to bytes
				const uint64_t permute_indexes = _pext_u64(0x0706050403020100, byte_mask);

				masks[i] = permute_indexes;
			}

			return masks;
		}();
		static const std::array<uint8_t, 16> out_idx_advance = [] {
			std::array<uint8_t, 16> lookup{};
			for (size_t i = 0; i < 16; ++i)
				lookup[i] = uint8_t(4 - pop_count(i));
			return lookup;
		}();

		inline_toggle static size_t pack_nonzero_sieve_chunks(const uint64_t* const sieve_data)
		{
			constexpr size_t trailing_chunks = (n_sieve_chunks % 4);
			constexpr size_t n_sieve_chunks_rounded = n_sieve_chunks - trailing_chunks;

			uint128_t indexes = _mm_setr_epi16(0, 1, 2, 3, 0, 0, 0, 0);
			const uint128_t inc = _mm_setr_epi16(4, 4, 4, 4, 0, 0, 0, 0);

			size_t out_idx = 0;

			// load one iteration ahead
			uint256_t data = _mm256_load_si256((uint256_t*)(&sieve_data[0]));

			for (chunk_idx_t i = 0; i < n_sieve_chunks_rounded; i += 4)
			{
				uint256_t mask_256 = _mm256_cmpeq_epi64(data, _mm256_setzero_si256()); // find the zero chunks
				const auto bit_mask = _mm256_movemask_pd(_mm256_castsi256_pd(mask_256));

				// convert element mask to shuffle/permute mask
				const uint128_t permute_mask = _mm_cvtsi64_si128(permute_masks[bit_mask]);

				// pack chunks
				const uint256_t packed_data = _mm256_permutevar8x32_epi32(data, _mm256_cvtepu8_epi32(permute_mask));

				// load one iteration ahead
				data = _mm256_load_si256((uint256_t*)(&sieve_data[i + 4]));

				// pack indexes
				const uint128_t keep_indexes = _mm_shuffle_epi8(indexes, permute_mask);

				// store packed chunks and indexes
				_mm256_storeu_si256((uint256_t*)(&sorted_chunks[0][out_idx]), packed_data);
				*(uint64_t*)(&chunk_indexes[0][out_idx]) = _mm_extract_epi64(keep_indexes, 0);

				out_idx += out_idx_advance[bit_mask]; // advance based on the number of elements we stored
				indexes = _mm_add_epi16(indexes, inc);
			}

			if constexpr (trailing_chunks >= 1)
			{
				const uint64_t chunk = sieve_data[n_sieve_chunks_rounded];
				sorted_chunks[0][out_idx] = chunk;
				chunk_indexes[0][out_idx] = n_sieve_chunks_rounded;
				out_idx += (chunk != 0);
			}
			if constexpr (trailing_chunks >= 2)
			{
				const uint64_t chunk = sieve_data[n_sieve_chunks_rounded + 1];
				sorted_chunks[0][out_idx] = chunk;
				chunk_indexes[0][out_idx] = n_sieve_chunks_rounded + 1;
				out_idx += (chunk != 0);
			}
			if constexpr (trailing_chunks == 3)
			{
				const uint64_t chunk = sieve_data[n_sieve_chunks_rounded + 2];
				sorted_chunks[0][out_idx] = chunk;
				chunk_indexes[0][out_idx] = n_sieve_chunks_rounded + 2;
				out_idx += (chunk != 0);
			}

			return out_idx; // the final index is the number of non-zero chunks
		}

		inline_toggle static void sort_chunks(const size_t n_nonzero_chunks)
		{
			using namespace detail;

			constexpr static uint32_t static_identity[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
			constexpr static uint8_t static_pc_shuf_lookup[32] = {
				0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
					0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
			constexpr static uint64_t static_nybble_mask[4] = {
				0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F };

			const uint64_t* const sieve_data = sorted_chunks[0].data();

			const size_t n_chunks_rounded = n_nonzero_chunks - (n_nonzero_chunks % 4);

			if (n_chunks_rounded >= 4)
			{
				alignas(32) uint32_t buffer[8]{};

				const uint256_t identity = _mm256_loadu_si256((uint256_t*)static_identity);
				const uint256_t seven_register = _mm256_set1_epi32(7);
				const uint256_t nybble_mask = _mm256_loadu_si256((uint256_t*)static_nybble_mask);
				const uint256_t pc_shuf_lookup = _mm256_loadu_si256((uint256_t*)static_pc_shuf_lookup);
				const uint256_t ymm_n_sieve_chunks = _mm256_set1_epi32(n_sieve_chunks);
				uint256_t pc_counts{};

				// run vector instructions one iteration ahead
				{
					uint256_t four_chunks = _mm256_loadu_si256((const uint256_t*)sieve_data);

					// sieve data -> nybbles
					const uint256_t nybbles_lo = _mm256_and_si256(four_chunks, nybble_mask);
					const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(four_chunks, 4), nybble_mask);
					// nybbles -> 8-bit pcs
					uint256_t popcounts = _mm256_add_epi8(_mm256_shuffle_epi8(pc_shuf_lookup, nybbles_lo),
														  _mm256_shuffle_epi8(pc_shuf_lookup, nybbles_hi));
					// 8-bit pcs -> 64-bit pcs
					popcounts = _mm256_sad_epu8(popcounts, _mm256_setzero_si256());

					// cap 4 popcounts to 7
					popcounts = _mm256_min_epu32(popcounts, seven_register); // using the 32-bit variant here is fine

					// calculate each chunk's offset into the first dimension of the array (popcount * n_sieve_chunks);
					uint256_t write_offsets = _mm256_mul_epu32(popcounts, ymm_n_sieve_chunks); // 32-bit multiply is okay here; half of the values are 0

					// duplicate each popcount so we can use four 4x64-bit permutes to produce eight copies of each popcount
					popcounts = _mm256_castps_si256(_mm256_moveldup_ps(_mm256_castsi256_ps(popcounts)));

					// chunk 1
					const uint256_t bcast_a = _mm256_permute4x64_epi64(popcounts, 0b00000000);
					const uint256_t shuffled_a = _mm256_permutevar8x32_epi32(pc_counts, bcast_a); // extract the existing count
					uint256_t mask = _mm256_cmpeq_epi32(bcast_a, identity); // generate a ones-mask at index [pc]
					pc_counts = _mm256_sub_epi32(pc_counts, mask); // subtract (-1) to increment the existing count

					// chunk 2
					const uint256_t bcast_b = _mm256_permute4x64_epi64(popcounts, 0b01010101);
					const uint256_t shuffled_b = _mm256_permutevar8x32_epi32(pc_counts, bcast_b);
					mask = _mm256_cmpeq_epi32(bcast_b, identity);
					uint256_t minor_offsets = _mm256_blend_epi32(shuffled_a, shuffled_b, 0b00001100);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					// chunk 3
					const uint256_t bcast_c = _mm256_permute4x64_epi64(popcounts, 0b10101010);
					const uint256_t shuffled_c = _mm256_permutevar8x32_epi32(pc_counts, bcast_c);
					mask = _mm256_cmpeq_epi32(bcast_c, identity);
					minor_offsets = _mm256_blend_epi32(minor_offsets, shuffled_c, 0b00110000);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					// chunk 4
					const uint256_t bcast_d = _mm256_permute4x64_epi64(popcounts, 0b11111111);
					const uint256_t shuffled_d = _mm256_permutevar8x32_epi32(pc_counts, bcast_d);
					mask = _mm256_cmpeq_epi32(bcast_d, identity);
					minor_offsets = _mm256_blend_epi32(minor_offsets, shuffled_d, 0b11000000);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					// add the second dimension of offsets
					write_offsets = _mm256_add_epi32(write_offsets, minor_offsets);

					// store results on the stack
					_mm256_storeu_si256((uint256_t*)buffer, write_offsets);
				}

				// load ahead
				uint256_t four_chunks = _mm256_loadu_si256((const uint256_t*)(sieve_data + 4));

				size_t i = 0;
				for (; i < n_chunks_rounded - 4; i += 4)
				{
					// run vector instructions one iteration ahead

					const uint256_t nybbles_lo = _mm256_and_si256(four_chunks, nybble_mask);
					const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(four_chunks, 4), nybble_mask);
					uint256_t popcounts = _mm256_add_epi8(_mm256_shuffle_epi8(pc_shuf_lookup, nybbles_lo),
														  _mm256_shuffle_epi8(pc_shuf_lookup, nybbles_hi));
					popcounts = _mm256_sad_epu8(popcounts, _mm256_setzero_si256());

					popcounts = _mm256_min_epu32(popcounts, seven_register);

					uint256_t write_offsets = _mm256_mul_epu32(popcounts, ymm_n_sieve_chunks);

					popcounts = _mm256_castps_si256(_mm256_moveldup_ps(_mm256_castsi256_ps(popcounts)));

					const uint256_t bcast_a = _mm256_permute4x64_epi64(popcounts, 0b00000000);
					const uint256_t shuffled_a = _mm256_permutevar8x32_epi32(pc_counts, bcast_a);
					uint256_t mask = _mm256_cmpeq_epi32(bcast_a, identity);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					const uint256_t bcast_b = _mm256_permute4x64_epi64(popcounts, 0b01010101);
					const uint256_t shuffled_b = _mm256_permutevar8x32_epi32(pc_counts, bcast_b);
					mask = _mm256_cmpeq_epi32(bcast_b, identity);
					uint256_t minor_offsets = _mm256_blend_epi32(shuffled_a, shuffled_b, 0b00001100);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					const uint256_t bcast_c = _mm256_permute4x64_epi64(popcounts, 0b10101010);
					const uint256_t shuffled_c = _mm256_permutevar8x32_epi32(pc_counts, bcast_c);
					mask = _mm256_cmpeq_epi32(bcast_c, identity);
					minor_offsets = _mm256_blend_epi32(minor_offsets, shuffled_c, 0b00110000);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					const uint256_t bcast_d = _mm256_permute4x64_epi64(popcounts, 0b11111111);
					const uint256_t shuffled_d = _mm256_permutevar8x32_epi32(pc_counts, bcast_d);
					mask = _mm256_cmpeq_epi32(bcast_d, identity);
					minor_offsets = _mm256_blend_epi32(minor_offsets, shuffled_d, 0b11000000);
					pc_counts = _mm256_sub_epi32(pc_counts, mask);

					write_offsets = _mm256_add_epi32(write_offsets, minor_offsets);


					// load two iterations ahead
					four_chunks = _mm256_loadu_si256((const uint256_t*)(sieve_data + i + 8));


					const uint32_t idx_0 = buffer[0];
					sorted_chunks.data()->_Elems[idx_0] = sieve_data[i + 0];
					chunk_indexes.data()->_Elems[idx_0] = chunk_indexes[0][i + 0];

					const uint32_t idx_1 = buffer[2];
					sorted_chunks.data()->_Elems[idx_1] = sieve_data[i + 1];
					chunk_indexes.data()->_Elems[idx_1] = chunk_indexes[0][i + 1];

					const uint32_t idx_2 = buffer[4];
					sorted_chunks.data()->_Elems[idx_2] = sieve_data[i + 2];
					chunk_indexes.data()->_Elems[idx_2] = chunk_indexes[0][i + 2];

					const uint32_t idx_3 = buffer[6];
					sorted_chunks.data()->_Elems[idx_3] = sieve_data[i + 3];
					chunk_indexes.data()->_Elems[idx_3] = chunk_indexes[0][i + 3];


					// store above results on the stack for the next iteration
					_mm256_storeu_si256((uint256_t*)buffer, write_offsets);
				}

				// cleanup step for the last unrolled iteration
				{
					const uint32_t idx_0 = buffer[0];
					sorted_chunks.data()->_Elems[idx_0] = sieve_data[i + 0];
					chunk_indexes.data()->_Elems[idx_0] = chunk_indexes[0][i + 0];

					const uint32_t idx_1 = buffer[2];
					sorted_chunks.data()->_Elems[idx_1] = sieve_data[i + 1];
					chunk_indexes.data()->_Elems[idx_1] = chunk_indexes[0][i + 1];

					const uint32_t idx_2 = buffer[4];
					sorted_chunks.data()->_Elems[idx_2] = sieve_data[i + 2];
					chunk_indexes.data()->_Elems[idx_2] = chunk_indexes[0][i + 2];

					const uint32_t idx_3 = buffer[6];
					sorted_chunks.data()->_Elems[idx_3] = sieve_data[i + 3];
					chunk_indexes.data()->_Elems[idx_3] = chunk_indexes[0][i + 3];
				}

				// save pc counts so far
				_mm256_storeu_si256((uint256_t*)buffer, pc_counts);
				for (size_t idx = 0; idx < 8; ++idx)
				{
					n_chunks_with_pc[idx] = chunk_count_t(buffer[idx]);
				}
			}

			// handle 0-3 remaining elements
			for (size_t i = n_chunks_rounded; i < n_nonzero_chunks; ++i)
			{
				const uint64_t chunk = sieve_data[i];

				size_t pc = (size_t)pop_count(chunk);
				pc = util::min(pc, 7ull);

				const chunk_count_t idx = n_chunks_with_pc[pc]++;

				sorted_chunks[pc][idx] = chunk;
				chunk_indexes[pc][idx] = chunk_indexes[0][i];
			}
		}

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
		inline_toggle void extract_candidates_with_popcount(uint64_t*& candidates)
		{
			const size_t n_chunks = n_chunks_with_pc[popcount];
			size_t i = 0;

			if constexpr (popcount <= 2)
			{
				const size_t n_chunks_rounded = n_chunks - (n_chunks % 4);
			#pragma clang loop vectorize(disable)
				for (; i < n_chunks_rounded; i += 4)
				{
					uint64_t chunk = sorted_chunks[popcount][i];
					uint64_t index = chunk_indexes[popcount][i] * 64ull;
					extract_candidates<popcount>(chunk, index, candidates);
					chunk = sorted_chunks[popcount][i + 1];
					index = chunk_indexes[popcount][i + 1] * 64ull;
					extract_candidates<popcount>(chunk, index, candidates);
					chunk = sorted_chunks[popcount][i + 2];
					index = chunk_indexes[popcount][i + 2] * 64ull;
					extract_candidates<popcount>(chunk, index, candidates);
					chunk = sorted_chunks[popcount][i + 3];
					index = chunk_indexes[popcount][i + 3] * 64ull;
					extract_candidates<popcount>(chunk, index, candidates);
				}
			}

			for (; i < n_chunks; ++i)
			{
				uint64_t chunk = sorted_chunks[popcount][i];
				uint64_t index = chunk_indexes[popcount][i] * 64ull;

				// generate instructions to extract exactly n candidates
				extract_candidates<popcount>(chunk, index, candidates);
			}
		}

		inline_toggle static void extract_bit_indexes(uint64_t*& candidates)
		{
			extract_candidates_with_popcount<1>(candidates);
			extract_candidates_with_popcount<2>(candidates);
			extract_candidates_with_popcount<3>(candidates);
			extract_candidates_with_popcount<4>(candidates);
			extract_candidates_with_popcount<5>(candidates);
			extract_candidates_with_popcount<6>(candidates);

			const size_t n_chunks_with_seven_or_more_bits = n_chunks_with_pc[7];
			for (size_t j = 0; j < n_chunks_with_seven_or_more_bits; ++j)
			{
				uint64_t chunk = sorted_chunks[7][j];
				const uint64_t index = chunk_indexes[7][j] * 64ull;

				*candidates++ = index + _tzcnt_u64(chunk); // store the chunk's index plus the next bit's index
				chunk = _blsr_u64(chunk); // reset the bit we just read
				*candidates++ = index + _tzcnt_u64(chunk); // bit 2
				chunk = _blsr_u64(chunk);
				*candidates++ = index + _tzcnt_u64(chunk); // 3
				chunk = _blsr_u64(chunk);
				*candidates++ = index + _tzcnt_u64(chunk); // 4
				chunk = _blsr_u64(chunk);
				*candidates++ = index + _tzcnt_u64(chunk); // 5
				chunk = _blsr_u64(chunk);
				*candidates++ = index + _tzcnt_u64(chunk); // 6
				chunk = _blsr_u64(chunk);

				do // bits 7+
				{
					*candidates++ = index + _tzcnt_u64(chunk);
					chunk = _blsr_u64(chunk);
				} while (chunk != 0);

			} // end for (each chunk with 7 or more bits)
		}

		inline_toggle static void convert_indexes_to_bitstrings(uint64_t* const candidates_start,
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
	inline_toggle static uint64_t* gather_sieve_results(uint64_t* candidates,
														const sieve_container& sieve,
														const uint64_t number)
	{
		using namespace detail;

		// autovectorizes to three large writes
		for (chunk_count_t& pc : n_chunks_with_pc)
			pc = 0;

		// pack non-zero sieve chunks before bitmap decoding
		const size_t n_nonzero_chunks = pack_nonzero_sieve_chunks((const uint64_t* const)sieve.data());

		// prepare arrays of 64-bit chunks, where the chunks in array n contain exactly n set bits (n candidates)
		sort_chunks(n_nonzero_chunks);

		uint64_t* const candidates_start = candidates;
		extract_bit_indexes(candidates);

		convert_indexes_to_bitstrings(candidates_start, candidates, number);

		return candidates;
	}

	namespace detail
	{
		void verify_sieve_offset_cache(const uint64_t start);
	}
}
