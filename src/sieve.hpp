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

			// We now have the distance to the next odd multiple of p.
			// Divide by 2 to store the *index* of the next odd multiple of p.
			sieve_offsets_cache[i] = n / 2;
		}
	}



	void partial_sieve(sieve_container& sieve
					   count_passes(, size_t& ps15))
	{
		// Start with the first prime not in the static sieve
		const sieve_prime_t* prime_ptr = small_primes_lookup.data() + static_sieve_primes.size() + 1;
		sieve_offset_t* offset_cache_ptr = sieve_offsets_cache.data() + static_sieve_primes.size() + 1;

		constexpr size_t sieve_end = sieve.size();

		// Sieve primes by strides of 15*p:
		// 
		// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
		//    x  x     x        x  x        x     x  x   <-- values to mark composite/false
		// x        x     x  x        x  x     x         <-- values to ignore

		size_t next_p = *prime_ptr;
		size_t next_offset = *offset_cache_ptr;

		for (;;)
		{
			// Get the next prime, loading one iteration ahead
			const size_t p = next_p;
			next_p = *(prime_ptr + 1);

			// Get the position of the next odd multiple of p*15, loading one iteration ahead
			size_t j = next_offset;
			next_offset = *(offset_cache_ptr + 1);

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
			*offset_cache_ptr = sieve_offset_t((j + (15 * p)) - sieve_end);

			++prime_ptr;
			++offset_cache_ptr;

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
		constexpr size_t block_size = ((sieve.size() / 4) / 8) * 8;
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
