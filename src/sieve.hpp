#pragma once

#include "util/types.hpp"

namespace mbp
{
	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(),
														 static_sieve_primes.end(), size_t(1), std::multiplies());

	using sieve_t = uint8_t;
	using sieve_container = std::array<sieve_t, static_sieve_size>;
	consteval sieve_container generate_static_sieve()
	{
		// sieve_container sieve(static_sieve_size, true);
		sieve_container sieve{}; std::fill(begin(sieve), end(sieve), true);

		// for each prime, mark off all multiples
		for (const auto p : static_sieve_primes)
			for (size_t i = 0; i < sieve.size(); i += p)
				sieve[i] = false;

		return sieve;
	}
	static alignas(64) constexpr sieve_container static_sieve = generate_static_sieve();

	// When sieving in steps of (some factor N) * (some prime P), we may have skipped a few multiples of P on
	// on either end of the sieve. We handle these with a set of branchless writes on each end of the sieve.
	// This overhead has a constant cost, which we can use to estimate the threshold at which it is no longer
	// an optimization to pass (larger) small primes through the stepped loop. This threshold is calculated at
	// compile time, based on requiring primes to make at least X writes, where X is a tuneable knob.
	consteval sieve_prime_t get_threshold(const size_t X)
	{
		for (size_t i = 1; i < small_primes_lookup.size(); ++i)
			if (small_primes_lookup[i] > (static_sieve_size / X))
				return small_primes_lookup[i - 1];
		// compile-time error if we don't find a valid answer
	}
	constexpr sieve_prime_t last_prime_for_stepping_by_fifteen = get_threshold((5 + 1) * 15);

	using sieve_offset_t = util::narrowest_uint_for_val<static_sieve_size>;
	std::vector<sieve_offset_t> sieve_offsets_cache(small_primes_lookup.size());

	void set_up_sieve_offsets_cache(const size_t start)
	{
		// Start with the first prime not in the static sieve.
		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			sieve_prime_t p = small_primes_lookup[i];

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
			sieve_prime_t n = p - (start % p);

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



	void partial_sieve_by_15(sieve_container& sieve,
							 const sieve_prime_t*& prime_ptr,
							 sieve_prime_t*& offset_cache_ptr)
	{
		sieve_t* sieve_begin = sieve.data();
		const sieve_t* const sieve_end = sieve_begin + sieve.size();

		size_t next_p = *prime_ptr;
		sieve_t* next_offset = sieve_begin + *offset_cache_ptr;

		for (;;)
		{
			// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
			//    x  x     x        x  x        x     x  x   <-- values to mark composite/false
			// x        x     x  x        x  x     x         <-- values to ignore

			// Get the next prime, loading one iteration ahead
			const size_t p = next_p;
			next_p = *(prime_ptr + 1);

			// Get the position of the next odd multiple of p*15, loading one iteration ahead
			sieve_t* j = next_offset;
			next_offset = sieve_begin + *(offset_cache_ptr + 1);

			// j is aligned to a multiple of 15, so we have 0-8 multiples of p *behind* j that must be crossed off.
			// Use eight branchless writes to j-(N * p) if they exist, falling back to a write to j in each case.
			*(((j - (14 * p)) >= sieve_begin) ? (j - (14 * p)) : j) = false;
			*(((j - (13 * p)) >= sieve_begin) ? (j - (13 * p)) : j) = false;
			*(((j - (11 * p)) >= sieve_begin) ? (j - (11 * p)) : j) = false;
			*(((j - (8 * p)) >= sieve_begin) ? (j - (8 * p)) : j) = false;
			*(((j - (7 * p)) >= sieve_begin) ? (j - (7 * p)) : j) = false;
			*(((j - (4 * p)) >= sieve_begin) ? (j - (4 * p)) : j) = false;
			*(((j - (2 * p)) >= sieve_begin) ? (j - (2 * p)) : j) = false;
			*(((j - p) >= sieve_begin) ? (j - p) : j) = false;

			// Stop marking 15p early, then handle cleanup after the loop.
			const sieve_t* const padded_end = sieve_end - (15 * p);

			// Each iteration, j points to an (already marked) multiple of p*15.
			// Mark off specific, implicitly odd, multiples of p.
			do
			{
				// skip 0
				*(j + p) = false;
				*(j + (2 * p)) = false;
				// skip three
				*(j + (4 * p)) = false;
				// skip five
				// skip six
				*(j + (7 * p)) = false;
				*(j + (8 * p)) = false;
				// skip nine
				// skip 10
				*(j + (11 * p)) = false;
				// skip 12
				*(j + (13 * p)) = false;
				*(j + (14 * p)) = false;

				j += (15 * p);
			} while (j < padded_end);

			// Same as above, perform any remaining writes
			*(((j + p) < sieve_end) ? (j + p) : j) = false;
			*(((j + (2 * p)) < sieve_end) ? (j + (2 * p)) : j) = false;
			*(((j + (4 * p)) < sieve_end) ? (j + (4 * p)) : j) = false;
			*(((j + (7 * p)) < sieve_end) ? (j + (7 * p)) : j) = false;
			*(((j + (8 * p)) < sieve_end) ? (j + (8 * p)) : j) = false;
			*(((j + (11 * p)) < sieve_end) ? (j + (11 * p)) : j) = false;
			*(((j + (13 * p)) < sieve_end) ? (j + (13 * p)) : j) = false;
			*(((j + (14 * p)) < sieve_end) ? (j + (14 * p)) : j) = false;

			// Calculate the offset and update the cache for the next sieving
			*offset_cache_ptr = sieve_prime_t((j + (15 * p)) - sieve_end);

			++prime_ptr;
			++offset_cache_ptr;

			if (p == last_prime_for_stepping_by_fifteen) break;
		}
	}

	void partial_sieve_by_3(sieve_container& sieve,
							const sieve_prime_t*& prime_ptr,
							sieve_prime_t*& offset_cache_ptr)
	{
		sieve_t* sieve_begin = sieve.data();
		const sieve_t* const sieve_end = sieve_begin + sieve.size();
		const sieve_prime_t* const small_primes_end = small_primes_lookup.data() + small_primes_lookup.size();

		size_t next_p = *prime_ptr;
		sieve_t* next_offset = sieve_begin + *offset_cache_ptr;

		for (; prime_ptr < small_primes_end; ++prime_ptr, ++offset_cache_ptr)
		{
			const size_t p = next_p;
			next_p = *(prime_ptr + 1);

			sieve_t* j = next_offset;
			next_offset = sieve_begin + *(offset_cache_ptr + 1);

			// Same as above, we are stepping by multiples of p. Use two branchless writes
			// before and after the loop to handle any multiples near the start and end of the sieve.
			*(((j - (2 * p)) >= sieve_begin) ? (j - (2 * p)) : j) = false;
			*(((j - p) >= sieve_begin) ? (j - p) : j) = false;

			const sieve_t* const padded_end = sieve_end - (3 * p);

			// Same as above, j points to an (already marked) multiple of p.
			// Only mark the necessary multiples of p.
			do
			{
				*(j + p) = false;
				*(j + (2 * p)) = false;
				j += (3 * p);
			} while (j < padded_end);

			*(((j + p) < sieve_end) ? (j + p) : j) = false;
			*(((j + (2 * p)) < sieve_end) ? (j + (2 * p)) : j) = false;

			*offset_cache_ptr = sieve_prime_t((j + (3 * p)) - sieve_end);
		}
	}

	void partial_sieve(sieve_container& sieve
					   count_passes(, size_t& ps15, size_t& ps3))
	{
		// Start with the first prime not in the static sieve
		const sieve_prime_t* prime_ptr = small_primes_lookup.data() + static_sieve_primes.size() + 1;
		sieve_prime_t* offset_cache_ptr = sieve_offsets_cache.data() + static_sieve_primes.size() + 1;

		// Sieve small primes by strides of 15*p
		partial_sieve_by_15(sieve, prime_ptr, offset_cache_ptr);
		count_passes(ps15 += util::vector_count_ones(sieve.data(), sieve.size()));

		// Sieve larger primes by smaller strides of 3*p
		partial_sieve_by_3(sieve, prime_ptr, offset_cache_ptr);
		count_passes(ps3 += util::vector_count_ones(sieve.data(), sieve.size()));
	}

	namespace detail
	{
		__forceinline void vectorized_gather_step(const sieve_t*& ptr, uint256_t& ymm,
												  size_t*& candidates, const size_t number,
												  const sieve_t* const sieve_begin)
		{
			ymm = _mm256_cmpeq_epi8(ymm, _mm256_set1_epi8(1));
			const uint32_t mask = (uint32_t)_mm256_movemask_epi8(ymm);
			const uint32_t trailing_zeroes = _tzcnt_u32(mask);
			const uint32_t mask_is_non_zero = !!mask; // map 0 -> 0; non-zero -> 1

			ptr += trailing_zeroes + mask_is_non_zero; // advance ptr and number by (tzc + 1) or 32, whichever is smaller
			ymm = _mm256_loadu_si256((uint256_t*)ptr); // load next block early

			*candidates = number + (2 * (ptr - sieve_begin)) - 2; // always store to avoid a branch			
			candidates += mask_is_non_zero; // keep the candidate (increment the pointer) unless mask is 0
		}

		__forceinline void scalar_gather(const sieve_t*& ptr, const sieve_t* const block_end,
										 size_t*& candidates, const size_t number,
										 const sieve_t* const sieve_begin)
		{
			for (; ptr < block_end; ++ptr)
			{
				*candidates = number + 2 * (ptr - sieve_begin);
				candidates += *ptr;
			}
		}
	}

	template<size_t n_bytes>
	size_t* gather_sieve_results_vectorized(size_t* candidates,
											const sieve_t* sieve_ptr,
											const size_t number)
	{
		static_assert(sizeof(sieve_t) == 1);

		constexpr size_t block_size_in_bytes = n_bytes / 4;

		const sieve_t* block_0 = sieve_ptr;
		const sieve_t* block_1 = sieve_ptr + (block_size_in_bytes * 1);
		const sieve_t* block_2 = sieve_ptr + (block_size_in_bytes * 2);
		const sieve_t* block_3 = sieve_ptr + (block_size_in_bytes * 3);

		const sieve_t* block_0_end = block_1;
		const sieve_t* block_1_end = block_2;
		const sieve_t* block_2_end = block_3;
		const sieve_t* block_3_end = sieve_ptr + n_bytes;

		// load ahead
		uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)block_0);
		uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)block_1);
		uint256_t ymm2 = _mm256_loadu_si256((uint256_t*)block_2);
		uint256_t ymm3 = _mm256_loadu_si256((uint256_t*)block_3);



		// Gather from 4 locations in parallel

		size_t safe_reads = block_size_in_bytes / 32;

		do
		{
			do
			{
				detail::vectorized_gather_step(block_0, ymm0, candidates, number, sieve_ptr);
				detail::vectorized_gather_step(block_1, ymm1, candidates, number, sieve_ptr);
				detail::vectorized_gather_step(block_2, ymm2, candidates, number, sieve_ptr);
				detail::vectorized_gather_step(block_3, ymm3, candidates, number, sieve_ptr);
			} while (--safe_reads != 0);

			safe_reads = util::min(block_0_end - block_0,
								   block_1_end - block_1,
								   block_2_end - block_2,
								   block_3_end - block_3) / 32;
		} while (safe_reads != 0);

		// At least one ptr is within 32 bytes of the end of its block.
		// If it is one of the first three blocks, swap with block 3 for cleanup

		if (block_0_end - block_0 < 32)
		{
			std::swap(block_0, block_3);
			std::swap(block_0_end, block_3_end);
			std::swap(ymm0, ymm3);
		}
		else if (block_1_end - block_1 < 32)
		{
			std::swap(block_1, block_3);
			std::swap(block_1_end, block_3_end);
			std::swap(ymm1, ymm3);
		}
		else if (block_2_end - block_2 < 32)
		{
			std::swap(block_2, block_3);
			std::swap(block_2_end, block_3_end);
			std::swap(ymm2, ymm3);
		}

		// in all cases, finish whatever block 3 is now
		detail::scalar_gather(block_3, block_3_end, candidates, number, sieve_ptr);



		// Gather from 3 locations in parallel

		// recalculate
		safe_reads = util::min(block_0_end - block_0,
							   block_1_end - block_1,
							   block_2_end - block_2) / 32;

		while (safe_reads != 0)
		{
			do
			{
				detail::vectorized_gather_step(block_0, ymm0, candidates, number, sieve_ptr);
				detail::vectorized_gather_step(block_1, ymm1, candidates, number, sieve_ptr);
				detail::vectorized_gather_step(block_2, ymm2, candidates, number, sieve_ptr);
			} while (--safe_reads != 0);

			safe_reads = util::min(block_0_end - block_0,
								   block_1_end - block_1,
								   block_2_end - block_2) / 32;
		}

		if (block_0_end - block_0 < 32)
		{
			std::swap(block_0, block_2);
			std::swap(block_0_end, block_2_end);
			std::swap(ymm0, ymm2);
		}
		else if (block_1_end - block_1 < 32)
		{
			std::swap(block_1, block_2);
			std::swap(block_1_end, block_2_end);
			std::swap(ymm1, ymm2);
		}

		detail::scalar_gather(block_2, block_2_end, candidates, number, sieve_ptr);



		// Gather from 2 locations in parallel

		safe_reads = util::min(block_0_end - block_0,
							   block_1_end - block_1) / 32;

		while (safe_reads != 0)
		{
			do
			{
				detail::vectorized_gather_step(block_0, ymm0, candidates, number, sieve_ptr);
				detail::vectorized_gather_step(block_1, ymm1, candidates, number, sieve_ptr);
			} while (--safe_reads != 0);

			safe_reads = util::min(block_0_end - block_0,
								   block_1_end - block_1) / 32;
		}

		if (block_0_end - block_0 < 32)
		{
			std::swap(block_0, block_1);
			std::swap(block_0_end, block_1_end);
			std::swap(ymm0, ymm1);
		}

		detail::scalar_gather(block_1, block_1_end, candidates, number, sieve_ptr);



		// Gather from remaining pointer

		safe_reads = (block_0_end - block_0) / 32;

		while (safe_reads != 0)
		{
			do
			{
				detail::vectorized_gather_step(block_0, ymm0, candidates, number, sieve_ptr);
			} while (--safe_reads != 0);

			safe_reads = (block_0_end - block_0) / 32;
		}

		detail::scalar_gather(block_0, block_0_end, candidates, number, sieve_ptr);



		return candidates;
	}

}
