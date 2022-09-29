
// Hacky way to hide the contents of zmmintrin.h
// MSVC always includes it from immintrin.h
#ifndef _ZMMINTRIN_H_INCLUDED
#define _ZMMINTRIN_H_INCLUDED
#endif

#include <bitset>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#define VCL_NAMESPACE vcl
#define MAX_VECTOR_SIZE 256
#include "../lib/vcl/vectorclass.h"
#include "config.hpp"
#include "io/io.hpp"
#include "math/franken_mpir.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/sandbox.hpp"
#include "util/simd.hpp"
#include "util/types.hpp"
#include "util/utility.hpp"


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
	static alignas(sieve_alignment) constexpr sieve_container static_sieve = generate_static_sieve();

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
	constexpr sieve_prime_t last_prime_for_stepping_by_threes = get_threshold(16);

	using sieve_offset_t = narrowest_uint_for_val<static_sieve_size>;
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
			else if (p <= last_prime_for_stepping_by_threes)
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

	void partial_sieve(sieve_container& sieve)
	{
		// *2 because we always take at least one step, and still need step*p padding
		static_assert(last_prime_for_stepping_by_fifteen < static_sieve_size / (15 * 2));
		static_assert(last_prime_for_stepping_by_fifteen < last_prime_for_stepping_by_threes);
		static_assert(last_prime_for_stepping_by_threes < static_sieve_size / (3 * 2));

		sieve_t* sieve_begin = sieve.data();
		const sieve_t* const sieve_end = sieve_begin + sieve.size();

		// Start with the first prime not in the static sieve
		const sieve_prime_t* prime_ptr = small_primes_lookup.data() + static_sieve_primes.size() + 1;
		sieve_prime_t* offset_cache_ptr = sieve_offsets_cache.data() + static_sieve_primes.size() + 1;

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

		// Sieve larger primes by smaller strides, using 3*p instead of 15*p

		for (;;)
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

			++prime_ptr;
			++offset_cache_ptr;

			if (p == last_prime_for_stepping_by_threes) break;
		}

		// Sieve the remaining primes, which are likely too large relative to the sieve size
		// to benefit from stepping by more than p.

		for (; prime_ptr < small_primes_lookup.data() + small_primes_lookup.size();
			 ++prime_ptr, ++offset_cache_ptr)
		{
			const size_t p = next_p;
			next_p = *(prime_ptr + 1);

			sieve_t* j = next_offset;
			next_offset = sieve_begin + *(offset_cache_ptr + 1);

			do
			{
				*j = false;
				j += p;
			} while (j < sieve_end);

			*offset_cache_ptr = sieve_prime_t(j - sieve_end);

			/*
			sieve_t* const j_start = j;

		#define try_next_write(n) \
			j = (j_start + p*n < sieve_end) ? (j_start + p*n) : j; \
			*j = false;
			// end define

			// 6 hardcoded writes
			*(j + p * 0) = false;
			*(j + p * 1) = false;
			*(j + p * 2) = false;
			*(j + p * 3) = false;
			*(j + p * 4) = false; // 5
			*(j + p * 5) = false;

			// Set j to the last written byte
			j += p * 5;

			// Use cmovs to make any further writes required
			try_next_write(6);
			try_next_write(7);
			try_next_write(8);
			try_next_write(9); // 10
			try_next_write(10);
			try_next_write(11);
			try_next_write(12);
			try_next_write(13);
			try_next_write(14); // 15
			try_next_write(15); // 16

			// Instead of stepping by p until we pass the end of the sieve, calculate at
			// compile time the max and min possible number of writes. Hardcode the min
			// possible writes, and use cmovs to make up the difference.
			constexpr size_t fewest_possible_writes_required = sieve_container().size() / small_primes_lookup.back();
			constexpr size_t most_possible_writes_required = sieve_container().size() / last_prime_for_stepping_by_threes;
			static_assert(fewest_possible_writes_required == 6); // Intellisense false-positive
			static_assert(most_possible_writes_required == 16); // Intellisense false-positive

			*offset_cache_ptr = sieve_prime_t((j + p) - sieve_end);
			*/
		}
	}

	tests_are_inlined size_t* gather_sieve_results(size_t* candidates,
												   const sieve_t* sieve_ptr, const sieve_t* const sieve_end,
												   size_t number)
	{
		static_assert(static_sieve_size % 15 == 0);

		// By scanning the sieve in steps of 3*5 == 15, we know every 3rd and every 5th value is false.
		// We can save work by never checking those ones:
		// 
		// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
		//    x  x     x        x  x        x     x  x   <-- values to check
		// x        x     x  x        x  x     x         <-- values to ignore
		//
		// In this case, we scan the sieve in steps of 3*5*7 == 105, which allows us to skip a few more
		// always-composite values, and only read the remaining 48 / 105.

		// "Offsets" stores the indexes that we must check, working from sieve_ptr.
		// sieve_ptr is aligned on a multiple of 3*5*7 from the beginning of the sieve,
		// which itself has stricter alignment on the number line.
		constexpr std::array<size_t, 48> offsets = []() consteval {
			std::array<size_t, 48> offsets{};
			size_t idx = 0;
			// Find the offsets that are not divisible by 3, 5, or 7
			for (size_t n = 1; n <= (3ull * 5 * 7); ++n)
				if (n % 3 != 0 && n % 5 != 0 && n % 7 != 0)
					offsets[idx++] = n;

			// The offset must be multiplied by 2 when calculating the bitstring, because the sieve only represents odd values.
			// There is a performance gotcha here: when this offset reaches 128, almost half of the encoded instructions below
			// become slightly larger, and no longer fit in 16 bytes. We can fix this by advancing both "number" and
			// "sieve_ptr" halfway through the loop, and then continuing with new, smaller offsets. We make these offsets smaller
			// by subtracting offset[23], the midpoint, from each of the following 24 offsets.
			for (size_t n = 24; n < offsets.size(); ++n)
				offsets[n] -= offsets[23];

			return offsets;
		}();

		for (; sieve_ptr < sieve_end;
			 sieve_ptr += (3ull * 5 * 7) - offsets[23],
			 number += (3ull * 5 * 7 * 2) - (offsets[23] * 2))
		{
			auto it = offsets.cbegin();

			*candidates = number + *it * 2; // unconditionally store number + offset*2
			candidates += *(sieve_ptr + *it++); // conditionally increment
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // write #5

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 10

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 15

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 20

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);

			number += *it * 2;
			sieve_ptr += *it++;
			*candidates = number;
			candidates += *(sieve_ptr); // 24 is a special case - see above

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 25

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 30

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 35

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 40

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 45

			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++);
			*candidates = number + *it * 2;
			candidates += *(sieve_ptr + *it++); // 48
		}

		return candidates;
	}



	tests_are_inlined size_t* prime_popcount_test(size_t* input,
												  const size_t* const candidates_end)
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();

		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) & 0b11);

		// Alternate between two output ptrs to split the workload
		// into two dependency chains.
		size_t* output_a = input;
		size_t* output_b = input + 1;

		for (; input < candidates_end_rounded; input += 4)
		{
			const size_t c0 = *input;
			const size_t pc_c0 = pop_count(c0);
			*output_a = c0;
			if (tiny_primes_lookup & (1ull << pc_c0)) output_a += 2;

			const size_t c1 = *(input + 1);
			const size_t pc_c1 = pop_count(c1);
			*output_b = c1;
			if (tiny_primes_lookup & (1ull << pc_c1)) output_b += 2;

			const size_t c2 = *(input + 2);
			const size_t pc_c2 = pop_count(c2);
			*output_a = c2;
			if (tiny_primes_lookup & (1ull << pc_c2)) output_a += 2;

			const size_t c3 = *(input + 3);
			const size_t pc_c3 = pop_count(c3);
			*output_b = c3;
			if (tiny_primes_lookup & (1ull << pc_c3)) output_b += 2;
		}

		// We've effectively created two interleaved arrays in one.
		// When the output ptrs end up nonadjacent, move data from a->b or b->a
		// until our data is contiguous again.

		// Scenario 1
		while (output_a > output_b + 1)
		{
			output_a -= 2;
			*output_b = *output_a;
			output_b += 2;
		}

		// Scenario 2
		while (output_a < output_b - 1)
		{
			output_b -= 2;
			*output_a = *output_b;
			output_a += 2;
		}

		// keep the smaller of the two
		size_t* output = (output_a < output_b) ? output_a : output_b;

		// handle last few elements
		for (; input < candidates_end; ++input)
		{
			const size_t c0 = *input;
			const size_t pc_c0 = pop_count(c0);
			*output = c0;
			if (tiny_primes_lookup & (1ull << pc_c0)) ++output;
		}

		return output;
	}

	tests_are_inlined size_t* gcd_test(size_t* input,
									   const size_t* const candidates_end)
	{
		constexpr size_t gcd_lookup = build_gcd_lookup();

		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) & 0b11);

		size_t* output = input;

		for (; input < candidates_end_rounded; input += 4)
		{
			const size_t n0 = *input;
			*output = n0;
			const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++output;

			const size_t n1 = *(input + 1);
			*output = n1;
			const size_t n1_pca = pop_count(n1 & 0xAAAAAAAAAAAAAAAA);
			const size_t n1_pcb = pop_count(n1 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n1_pca + 32 - n1_pcb))) ++output;

			const size_t n2 = *(input + 2);
			*output = n2;
			const size_t n2_pca = pop_count(n2 & 0xAAAAAAAAAAAAAAAA);
			const size_t n2_pcb = pop_count(n2 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n2_pca + 32 - n2_pcb))) ++output;

			const size_t n3 = *(input + 3);
			*output = n3;
			const size_t n3_pca = pop_count(n3 & 0xAAAAAAAAAAAAAAAA);
			const size_t n3_pcb = pop_count(n3 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n3_pca + 32 - n3_pcb))) ++output;
		}

		// handle last few elements
		for (; input < candidates_end; ++input)
		{
			const size_t n0 = *input;
			*output = n0;
			const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++output;
		}

		return output;
	}



	tests_are_inlined size_t* div_tests_with_four_rems(size_t* input,
													   const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<3, 5>::val;
		static_assert(bitmask == bitmask_for<5, 13>::val &&
					  bitmask == bitmask_for<8, 13>::val &&
					  bitmask == bitmask_for<4, 17>::val);
		static_assert(period_of<bitmask>::val == 4);

		const prime_lookup_t* const prime_factor_lookup_ptr = prime_factor_lookup.data();

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			const size_t pc_0 = pop_count(number & (bitmask << 0));
			size_t b3m5_rem = pc_0;
			size_t b8m13_rem = pc_0;
			size_t b5m13_rem = pc_0;
			size_t b4m17_rem = pc_0;

			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b3m5_rem += pc_1 * pow_mod<3, 1, 5>::rem;
			b5m13_rem += pc_1 * pow_mod<5, 1, 13>::rem;
			b8m13_rem += pc_1 * pow_mod<8, 1, 13>::rem;
			b4m17_rem += pc_1 * pow_mod<4, 1, 17>::rem;

			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b3m5_rem += pc_2 * pow_mod<3, 2, 5>::rem;
			b5m13_rem += pc_2 * pow_mod<5, 2, 13>::rem;
			b8m13_rem += pc_2 * pow_mod<8, 2, 13>::rem;
			b4m17_rem += pc_2 * pow_mod<4, 2, 17>::rem;

			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b3m5_rem += pc_3 * pow_mod<3, 3, 5>::rem;
			b5m13_rem += pc_3 * pow_mod<5, 3, 13>::rem;
			b8m13_rem += pc_3 * pow_mod<8, 3, 13>::rem;
			b4m17_rem += pc_3 * pow_mod<4, 3, 17>::rem;

			size_t merged_masks = 0;
			merged_masks |= (prime_factor_lookup_ptr[b3m5_rem] & (1ull << get_prime_index<5>::idx));
			merged_masks |= (prime_factor_lookup_ptr[b5m13_rem] & (1ull << get_prime_index<13>::idx));
			merged_masks |= (prime_factor_lookup_ptr[b8m13_rem] & (1ull << get_prime_index<13>::idx));
			merged_masks |= (prime_factor_lookup_ptr[b4m17_rem] & (1ull << get_prime_index<17>::idx));

			// Only advance the pointer if the number is still a candidate
			output += !merged_masks;
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_three_rems(size_t* input,
														const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<4, 7>::val;
		static_assert(bitmask == bitmask_for<3, 13>::val &&
					  bitmask == bitmask_for<9, 13>::val);
		static_assert(period_of<bitmask>::val == 3);

		const prime_lookup_t* const prime_factor_lookup_ptr = prime_factor_lookup.data();

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			const size_t pc_0 = pop_count(number & (bitmask << 0));
			size_t b4_m7_rem = pc_0;
			size_t b3_m13_rem = pc_0;
			size_t b9_m13_rem = pc_0;

			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b4_m7_rem += pc_1 * pow_mod<4, 1, 7>::rem;
			b3_m13_rem += pc_1 * pow_mod<3, 1, 13>::rem;
			b9_m13_rem += pc_1 * pow_mod<9, 1, 13>::rem;

			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b4_m7_rem += pc_2 * pow_mod<4, 2, 7>::rem;
			b3_m13_rem += pc_2 * pow_mod<3, 2, 13>::rem;
			b9_m13_rem += pc_2 * pow_mod<9, 2, 13>::rem;

			// Only advance the pointer if the number is still a candidate
			size_t merged_masks = 0;
			merged_masks |= (prime_factor_lookup_ptr[b4_m7_rem] & (1ull << get_prime_index<7>::idx));
			merged_masks |= (prime_factor_lookup_ptr[b3_m13_rem] & (1ull << get_prime_index<13>::idx));
			merged_masks |= (prime_factor_lookup_ptr[b9_m13_rem] & (1ull << get_prime_index<13>::idx));

			output += !merged_masks;
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_six_rems(size_t* input,
													  const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<3, 7>::val;
		static_assert(bitmask == bitmask_for<5, 7>::val &&
					  bitmask == bitmask_for<4, 13>::val &&
					  bitmask == bitmask_for<10, 13>::val);
		static_assert(period_of<bitmask>::val == 6);

		// constexpr uint128_t static_rems0 = { always all ones };
		constexpr static uint128_t static_rems1 = uint128_t{ .m128i_u32{
			pow_mod<3, 1, 7>::rem,
			pow_mod<5, 1, 7>::rem,
			pow_mod<4, 1, 13>::rem,
			pow_mod<10, 1, 13>::rem } };
		constexpr static uint128_t static_rems2 = uint128_t{ .m128i_u32{
			pow_mod<3, 2, 7>::rem,
			pow_mod<5, 2, 7>::rem,
			pow_mod<4, 2, 13>::rem,
			pow_mod<10, 2, 13>::rem } };
		constexpr static uint128_t static_rems3 = uint128_t{ .m128i_u32{
			pow_mod<3, 3, 7>::rem,
			pow_mod<5, 3, 7>::rem,
			pow_mod<4, 3, 13>::rem,
			pow_mod<10, 3, 13>::rem } };
		constexpr static uint128_t static_rems4 = uint128_t{ .m128i_u32{
			pow_mod<3, 4, 7>::rem,
			pow_mod<5, 4, 7>::rem,
			pow_mod<4, 4, 13>::rem,
			pow_mod<10, 4, 13>::rem } };
		constexpr static uint128_t static_rems5 = uint128_t{ .m128i_u32{
			pow_mod<3, 5, 7>::rem,
			pow_mod<5, 5, 7>::rem,
			pow_mod<4, 5, 13>::rem,
			pow_mod<10, 5, 13>::rem } };

		const auto xmm_rems1 = _mm_loadu_si128(&static_rems1);
		const auto xmm_rems2 = _mm_loadu_si128(&static_rems2);
		const auto xmm_rems3 = _mm_loadu_si128(&static_rems3);
		const auto xmm_rems4 = _mm_loadu_si128(&static_rems4);
		const auto xmm_rems5 = _mm_loadu_si128(&static_rems5);

		const prime_lookup_t* const prime_factor_lookup_ptr = prime_factor_lookup.data();

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			auto xmm0 = _mm_set1_epi32((int)pop_count(number & (bitmask << 0)));
			auto xmm1 = _mm_set1_epi32((int)pop_count(number & (bitmask << 1)));
			auto xmm2 = _mm_set1_epi32((int)pop_count(number & (bitmask << 2)));
			auto xmm3 = _mm_set1_epi32((int)pop_count(number & (bitmask << 3)));
			auto xmm4 = _mm_set1_epi32((int)pop_count(number & (bitmask << 4)));
			auto xmm5 = _mm_set1_epi32((int)pop_count(number & (bitmask << 5)));

			// multiply each popcount with Nth remainder of each test
			// (rems0 is always all ones)
			xmm1 = _mm_mullo_epi32(xmm1, xmm_rems1);
			xmm2 = _mm_mullo_epi32(xmm2, xmm_rems2);
			xmm3 = _mm_mullo_epi32(xmm3, xmm_rems3);
			xmm4 = _mm_mullo_epi32(xmm4, xmm_rems4);
			xmm5 = _mm_mullo_epi32(xmm5, xmm_rems5);

			// add pairs
			xmm0 = _mm_add_epi32(xmm0, xmm1);
			xmm2 = _mm_add_epi32(xmm2, xmm3);
			xmm4 = _mm_add_epi32(xmm4, xmm5);
			// final adds
			xmm0 = _mm_add_epi32(xmm0, xmm2);
			xmm0 = _mm_add_epi32(xmm0, xmm4);

			// Only advance the pointer if the number is still a candidate
			size_t merged_masks = 0;
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[0]] & (1ull << get_prime_index<7>::idx));
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[1]] & (1ull << get_prime_index<7>::idx));
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[2]] & (1ull << get_prime_index<13>::idx));
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[3]] & (1ull << get_prime_index<13>::idx));

			output += !merged_masks;
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_five_rems(size_t* input,
													   const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<3, 11>::val;  // base 3 % 11   (5 remainders)
		static_assert(bitmask == bitmask_for<4, 11>::val);  // base 4 % 11   (5 remainders)
		static_assert(bitmask == bitmask_for<5, 11>::val);  // base 5 % 11   (5 remainders)
		static_assert(bitmask == bitmask_for<9, 11>::val);  // base 9 % 11   (5 remainders)
		static_assert(period_of<bitmask>::val == 5);

		// constexpr uint128_t static_rems0 = { always all ones };
		constexpr static uint128_t static_rems1 = uint128_t{ .m128i_u32{
			pow_mod<3, 1, 11>::rem,
			pow_mod<4, 1, 11>::rem,
			pow_mod<5, 1, 11>::rem,
			pow_mod<9, 1, 11>::rem } };
		constexpr static uint128_t static_rems2 = uint128_t{ .m128i_u32{
			pow_mod<3, 2, 11>::rem,
			pow_mod<4, 2, 11>::rem,
			pow_mod<5, 2, 11>::rem,
			pow_mod<9, 2, 11>::rem } };
		constexpr static uint128_t static_rems3 = uint128_t{ .m128i_u32{
			pow_mod<3, 3, 11>::rem,
			pow_mod<4, 3, 11>::rem,
			pow_mod<5, 3, 11>::rem,
			pow_mod<9, 3, 11>::rem } };
		constexpr static uint128_t static_rems4 = uint128_t{ .m128i_u32{
			pow_mod<3, 4, 11>::rem,
			pow_mod<4, 4, 11>::rem,
			pow_mod<5, 4, 11>::rem,
			pow_mod<9, 4, 11>::rem } };

		const auto xmm_rems1 = _mm_loadu_si128(&static_rems1);
		const auto xmm_rems2 = _mm_loadu_si128(&static_rems2);
		const auto xmm_rems3 = _mm_loadu_si128(&static_rems3);
		const auto xmm_rems4 = _mm_loadu_si128(&static_rems4);

		const int* const prime_factor_lookup_ptr = (const int*)prime_factor_lookup.data();

		size_t* output = input;

		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;

			// always write
			*output = number;

			auto xmm0 = _mm_set1_epi32((int)pop_count(number & (bitmask << 0)));
			auto xmm1 = _mm_set1_epi32((int)pop_count(number & (bitmask << 1)));
			auto xmm2 = _mm_set1_epi32((int)pop_count(number & (bitmask << 2)));
			auto xmm3 = _mm_set1_epi32((int)pop_count(number & (bitmask << 3)));
			auto xmm4 = _mm_set1_epi32((int)pop_count(number & (bitmask << 4)));

			// multiply each popcount with Nth remainder of each test
			// (rems0 is always all ones)
			xmm1 = _mm_mullo_epi32(xmm1, xmm_rems1);
			xmm2 = _mm_mullo_epi32(xmm2, xmm_rems2);
			xmm3 = _mm_mullo_epi32(xmm3, xmm_rems3);
			xmm4 = _mm_mullo_epi32(xmm4, xmm_rems4);

			// add pairs
			xmm0 = _mm_add_epi32(xmm0, xmm1);
			xmm2 = _mm_add_epi32(xmm2, xmm3);
			// final adds
			xmm0 = _mm_add_epi32(xmm0, xmm2);
			xmm0 = _mm_add_epi32(xmm0, xmm4);

			size_t merged_lookups = prime_factor_lookup_ptr[xmm0.m128i_u32[0]];
			merged_lookups |= prime_factor_lookup_ptr[xmm0.m128i_u32[1]];
			merged_lookups |= prime_factor_lookup_ptr[xmm0.m128i_u32[2]];
			merged_lookups |= prime_factor_lookup_ptr[xmm0.m128i_u32[3]];

			constexpr size_t mask = (1ull << get_prime_index<11>::idx);

			// Only advance the pointer if each test's selected bit was 0
			output += !(merged_lookups & mask);

			// alternatively, use avx gather(), and(), and testz()
			//output += _mm_testz_si128(xmm0, prime_mask); // compute bitwise AND, return 1 if result is 0, else return 0
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_10_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// bases 6, 7, and 8 mod 11 (10 remainders)

		static_assert(sizeof(remainder_t) == 1);
		alignas(64) static constexpr remainder_t b6_rems[64] = {
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7 };
		alignas(64) static constexpr remainder_t b7_rems[64] = {
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2 };
		alignas(64) static constexpr remainder_t b8_rems[64] = {
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6 };

		const uint256_t b6_rems_lower = _mm256_loadu_si256((uint256_t*)&b6_rems[0]);
		const uint256_t b6_rems_upper = _mm256_loadu_si256((uint256_t*)&b6_rems[32]);

		const uint256_t b7_rems_lower = _mm256_loadu_si256((uint256_t*)&b7_rems[0]);
		const uint256_t b7_rems_upper = _mm256_loadu_si256((uint256_t*)&b7_rems[32]);

		const uint256_t b8_rems_lower = _mm256_loadu_si256((uint256_t*)&b8_rems[0]);
		const uint256_t b8_rems_upper = _mm256_loadu_si256((uint256_t*)&b8_rems[32]);

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			// Convert 64 bits to 64 bytes
			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));
			const uint256_t mask_upper = util::expand_bits_to_bytes(number >> 32);

			const uint256_t ymm0 = _mm256_and_si256(mask_lower, b6_rems_lower);
			const uint256_t ymm1 = _mm256_and_si256(mask_upper, b6_rems_upper);

			const uint256_t ymm2 = _mm256_and_si256(mask_lower, b7_rems_lower);
			const uint256_t ymm3 = _mm256_and_si256(mask_upper, b7_rems_upper);

			const uint256_t ymm4 = _mm256_and_si256(mask_lower, b8_rems_lower);
			const uint256_t ymm5 = _mm256_and_si256(mask_upper, b8_rems_upper);

			// Not all div tests can safely combine rems without moving to a larger type,
			// but our rems are small enough to avoid overflow.
			const auto b6_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm0, ymm1));
			const auto b7_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm2, ymm3));
			const auto b8_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm4, ymm5));

			size_t merged_lookups = prime_factor_lookup[b6_rem];
			merged_lookups |= prime_factor_lookup[b7_rem];
			merged_lookups |= prime_factor_lookup[b8_rem];

			// Only advance the pointer if the nth bit was 0 in all lookups
			if ((merged_lookups & (prime_lookup_t(1u) << get_prime_index<11>::idx)) == 0) ++output;
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_12_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		static_assert(sizeof(remainder_t) == 1);
		alignas(64) static constexpr remainder_t b6_rems[64] = {
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8 };
		alignas(64) static constexpr remainder_t b7_rems[64] = {
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5 };
		alignas(64) static constexpr remainder_t b11_rems[64] = {
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5 };

		const uint256_t b6_rems_lower = _mm256_loadu_si256((uint256_t*)&b6_rems[0]);
		const uint256_t b6_rems_upper = _mm256_loadu_si256((uint256_t*)&b6_rems[32]);

		const uint256_t b7_rems_lower = _mm256_loadu_si256((uint256_t*)&b7_rems[0]);
		const uint256_t b7_rems_upper = _mm256_loadu_si256((uint256_t*)&b7_rems[32]);

		const uint256_t b11_rems_lower = _mm256_loadu_si256((uint256_t*)&b11_rems[0]);
		const uint256_t b11_rems_upper = _mm256_loadu_si256((uint256_t*)&b11_rems[32]);

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			// Convert 64 bits to 64 bytes
			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));
			const uint256_t mask_upper = util::expand_bits_to_bytes(number >> 32);

			const uint256_t ymm0 = _mm256_and_si256(mask_lower, b6_rems_lower);
			const uint256_t ymm1 = _mm256_and_si256(mask_upper, b6_rems_upper);

			const uint256_t ymm2 = _mm256_and_si256(mask_lower, b7_rems_lower);
			const uint256_t ymm3 = _mm256_and_si256(mask_upper, b7_rems_upper);

			const uint256_t ymm4 = _mm256_and_si256(mask_lower, b11_rems_lower);
			const uint256_t ymm5 = _mm256_and_si256(mask_upper, b11_rems_upper);

			// Not all div tests can safely combine rems without moving to a larger type,
			// but our rems are small enough to avoid overflow.
			const auto b6_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm0, ymm1));
			const auto b7_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm2, ymm3));
			const auto b11_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm4, ymm5));

			prime_lookup_t merged_lookups = prime_factor_lookup[b6_rem];
			merged_lookups |= prime_factor_lookup[b7_rem];
			merged_lookups |= prime_factor_lookup[b11_rem];

			// Only advance the pointer if the nth bit was 0 in all lookups
			if ((merged_lookups & (prime_lookup_t(1u) << get_prime_index<13>::idx)) == 0) ++output;
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_16_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		constexpr static alignas(64) uint8_t static_rems[7 + 1][16] = {
			{ 1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6 }, // base 3 % 17
			{ 1, 5, 8, 6, 13, 14, 2, 10, 16, 12, 9, 11, 4, 3, 15, 7 }, // base 5 % 17
			{ 1, 6, 2, 12, 4, 7, 8, 14, 16, 11, 15, 5, 13, 10, 9, 3 }, // base 6 % 17
			{ 1, 7, 15, 3, 4, 11, 9, 12, 16, 10, 2, 14, 13, 6, 8, 5 }, // base 7 % 17
			{ 1, 10, 15, 14, 4, 6, 9, 5, 16, 7, 2, 3, 13, 11, 8, 12 }, // base 10 % 17
			{ 1, 11, 2, 5, 4, 10, 8, 3, 16, 6, 15, 12, 13, 7, 9, 14 }, // base 11 % 17
			{ 1, 12, 8, 11, 13, 3, 2, 7, 16, 5, 9, 6, 4, 14, 15, 10 }, // base 12 % 17
			{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } }; // padding

		const uint256_t rems_0_1 = _mm256_loadu_si256((uint256_t*)&static_rems[0]);
		const uint256_t rems_2_3 = _mm256_loadu_si256((uint256_t*)&static_rems[2]);
		const uint256_t rems_4_5 = _mm256_loadu_si256((uint256_t*)&static_rems[4]);
		const uint256_t rems_6_x = _mm256_loadu_si256((uint256_t*)&static_rems[6]);

		size_t* output = input;
		size_t number = *input;

		size_t upper_bits = number & upper_bits_mask;
		uint256_t upper_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number >> 32),
												 _mm256_set1_epi8(0x01));

		for (; input < candidates_end; )
		{
			// Upper 32 bits only change every 4B integers - re-calculate when stale
			if ((number & upper_bits_mask) != upper_bits)
			{
				upper_bits = number & upper_bits_mask;
				upper_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number >> 32),
											   _mm256_set1_epi8(0x01)); // convert 0xFF -> 0x01
			}

			uint256_t lower_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number & uint32_t(-1)),
													 _mm256_set1_epi8(0x01)); // convert 0xFF -> 0x01

			// always write
			*output = number;
			// load one iteration ahead
			++input;
			number = *input;

			// add the upper 32 bytes to the lower 32 bytes to get values in range 0-2
			uint256_t ymm0 = _mm256_add_epi8(lower_bytes, upper_bytes);

			// add the upper 16 bytes to the lower 16 bytes to get values in range 0-4, and leave the results in both lanes
			ymm0 = _mm256_add_epi8(ymm0, _mm256_permute2x128_si256(ymm0, ymm0, 1));

			// multiply 16+16 8-bit values by 32 remainders, storing partially summed results as 8+8 uint16_ts
			uint256_t ymm1 = _mm256_maddubs_epi16(ymm0, rems_0_1);
			uint256_t ymm2 = _mm256_maddubs_epi16(ymm0, rems_2_3);
			uint256_t ymm3 = _mm256_maddubs_epi16(ymm0, rems_4_5);
			uint256_t ymm4 = _mm256_maddubs_epi16(ymm0, rems_6_x);

			// pack four sparse registers into two
			ymm0 = _mm256_packus_epi16(ymm1, ymm2);
			ymm1 = _mm256_packus_epi16(ymm3, ymm4);

			// h-sum into 4 16-bit integers, 2 in each lane
			ymm0 = _mm256_sad_epu8(ymm0, _mm256_setzero_si256());
			ymm1 = _mm256_sad_epu8(ymm1, _mm256_setzero_si256());

			// extract high lanes
			const uint128_t xmm2 = _mm256_extracti128_si256(ymm0, 1);
			const uint128_t xmm3 = _mm256_extracti128_si256(ymm1, 1);

			prime_lookup_t merged_lookups = 0;
			merged_lookups |= prime_factor_lookup[ymm0.m256i_i64[0]];
			merged_lookups |= prime_factor_lookup[ymm0.m256i_i64[1]];
			merged_lookups |= prime_factor_lookup[ymm1.m256i_i64[0]];
			merged_lookups |= prime_factor_lookup[ymm1.m256i_i64[1]];
			merged_lookups |= prime_factor_lookup[xmm2.m128i_i64[0]];
			merged_lookups |= prime_factor_lookup[xmm2.m128i_i64[1]];
			merged_lookups |= prime_factor_lookup[xmm3.m128i_i64[0]];

			// Only advance the pointer if the nth bit was 0 in all lookups
			// Invert the bit we care about, move it to the 0th position, and mask to select it
			output += (~merged_lookups >> get_prime_index<17>::idx) & 0b1;
		}

		return output;
	}



	tests_are_inlined size_t* branchless_div_tests(size_t* const candidates_begin,
												   size_t* const candidates_end,
												   div_test::div_test_t* div_tests,
												   const size_t n_of_tests)
	{
		// See multibase_div_tests() for details.

		using namespace div_test;
		using namespace div_test::detail;
		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		size_t* shrinking_end = candidates_end;

		for (size_t i = 0; i < n_of_tests; ++i)
		{
			div_test_t& div_test = *(div_tests + i);

			const uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)&div_test.remainders[0]);
			const uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)&div_test.remainders[32]);

			size_t upper_bits = candidates_begin[0] & upper_bits_mask;
			const prime_lookup_t* pf_lookup_preinc = prime_factor_lookup.data() +
				util::vcl_hadd_x(_mm256_and_si256(util::expand_bits_to_bytes(upper_bits >> 32), ymm1));

			const size_t* input = candidates_begin;
			size_t* output = candidates_begin;

			size_t number = *input; // load one iteration ahead

			for (; input < shrinking_end; ++input)
			{
				*output = number; // always write

				// These only change every 4B integers - re-calculate when stale
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;
					pf_lookup_preinc = prime_factor_lookup.data() +
						util::vcl_hadd_x(_mm256_and_si256(util::expand_bits_to_bytes(upper_bits >> 32), ymm1));
				}

				const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));

				// load one iteration ahead
				number = *(input + 1);

				const uint256_t rems_lower = _mm256_and_si256(mask_lower, ymm0);
				const auto rem = util::vcl_hadd_x(rems_lower);

				// increment only if we haven't found a prime factor yet
				output += ((pf_lookup_preinc[rem] & (prime_lookup_t(1) << div_test.prime_idx)) == 0);
			}

			div_test.hits += uint32_t(input - output);

			// output is one past the last element; use it as the new end
			shrinking_end = output;
		}

		return shrinking_end;
	}

	tests_are_inlined size_t* multibase_div_tests(size_t* input,
												  const size_t* const candidates_end,
												  div_test::div_test_t* div_tests_start)
	{
		using namespace div_test;

		size_t* output = input;
		const auto* const div_tests_end = div_tests.data() + div_tests.size();

		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;

			// always write
			*output = number;

			static_assert(sizeof(remainder_t) == 1);

			// Convert each half of number to a 32-byte bitmask
			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));
			const uint256_t mask_upper = util::expand_bits_to_bytes(number >> 32);

			// Load 32+32 remainders into two 32-byte registers
			uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)&div_tests_start->remainders[0]);
			uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)&div_tests_start->remainders[32]);

			size_t is_candidate = 1;

			for (auto* div_test_ptr = div_tests_start; div_test_ptr < div_tests_end; ++div_test_ptr)
			{
				div_test_t& div_test = *div_test_ptr;

				// Use the byte-sized bits of the bitstring to select remainders
				const uint256_t rems_lower = _mm256_and_si256(mask_lower, ymm0);
				const uint256_t rems_upper = _mm256_and_si256(mask_upper, ymm1);

				// Load the next remainders, one iteration ahead
				ymm0 = _mm256_loadu_si256((uint256_t*)(((uint8_t*)&div_test.remainders) + sizeof(div_test_t) + 0));
				ymm1 = _mm256_loadu_si256((uint256_t*)(((uint8_t*)&div_test.remainders) + sizeof(div_test_t) + 32));

				// Calculate the horizontal sum of remainders. We have two vectors to h-sum,
				// but they can't be immediately added together without 8-bit overflow. We also
				// don't want to pay for two full h-sums. Solve this with a custom hadd: perform
				// the first hadd step on each vector, which extends their values from 8-bit to
				// 16-bit integers, then safely add the vectors together and continue the hadd.
				// This takes N+2 steps in total, instead of 2N.
				const size_t rem = util::vcl_hadd2_x(rems_upper, rems_lower);

				if (has_small_prime_factor(rem, div_test.prime_idx))
				{
					div_test.hits++;
					is_candidate = 0;
					break;
				}
			}

			// conditionally increment
			output += is_candidate;
		}

		return output;
	}

	void update_div_test_order()
	{
		using namespace div_test;

		for (size_t i = 1; i < div_tests.size(); ++i)
		{
			if (div_tests[i].hits > div_tests[i - 1].hits)
			{
				std::swap(div_tests[i], div_tests[i - 1]);

				// skip next comparison, so no test can move by more than one position per reorder
				++i;
			}
		}

		// Divide hit counts by 2 to create a weighted moving average
		for (auto& div_test : div_tests)
			div_test.hits >>= 1u;
	}

	void print_div_tests()
	{
	#if analyze_div_tests
		using namespace div_test;

		std::stringstream ss;

		ss << '\n';
		ss << "Prime factor lookup size: " << div_test::detail::prime_factor_lookup_size << '\n';
		ss << div_tests.size() << " div tests:\n";

		auto w = std::setw;
		for (const auto& dt : div_tests)
		{
			ss << "   base " << std::setfill(' ') << w(2) << size_t(dt.base) << " % " << w(3) << size_t(small_primes_lookup[dt.prime_idx]) << ": ";
			if (dt.hits == 0)
			{
				ss << "        -       ";
			}
			else
			{
				ss << w(9) << dt.hits << " hits  ";
			}

			ss << w(2) << size_t(dt.n_of_remainders) << " remainders: 1";
			for (size_t j = 1; j < dt.n_of_remainders; ++j)
			{
				if (j < 20)
				{
					ss << ' ' << w(3) << size_t(dt.remainders[j]);
				}
				else
				{
					ss << " ...";
					break;
				}
			}
			ss << '\n';
		}

		std::cout << ss.str();

	#endif
	}

	void run_div_test_analysis()
	{
	#if analyze_div_tests
		using namespace div_test;

		const auto div_test_pred = [](const auto& a, const auto& b) {
			return a.hits < b.hits;
		};

		if (std::is_sorted(div_tests.rbegin(), div_tests.rend(), div_test_pred))
		{
			std::cout << "Div tests have not changed hit count order\n";
			return;
		}

		const auto copy = div_tests;

		// sort descending
		std::sort(div_tests.rbegin(), div_tests.rend(), div_test_pred);

		size_t moved = 0;
		for (size_t i = 0; i < div_tests.size(); ++i)
			if (div_tests[i].base != copy[i].base || div_tests[i].prime_idx != copy[i].prime_idx)
				moved++;

		print_div_tests();

		static std::stringstream ss;
		ss << ' ' << moved;

		std::cout << moved << " div tests changed position\n";
		std::cout << '(' << ss.str() << ")\n\n";
	#endif
	}



	void find_multibase_primes()
	{
		gmp_randclass r{ gmp_randinit_mt };
		r.seed(mpir_ui{ 0xdeadbeef });

		size_t number = benchmark_mode ? bm_start : load_from_results();
		mpz_class mpz_number = 0ull; // it's a surprise tool that will help us later

		static alignas(sieve_alignment) sieve_container sieve {};

		// Round starting number down to the nearest odd multiple of the sieve sieze
		number -= sieve.size(); // n -= k
		number -= number % (2 * sieve.size()); // n -= n % 2k
		number += sieve.size(); // n += k

		set_up_sieve_offsets_cache(number);

		size_t next_div_test_reorder = number + div_test::reorder_interval;

		// 2x the expected number of candidates from the sieve passes
		constexpr size_t candidates_capacity = [] {
			double cleared = 0.0;
			for (size_t i = 1; i < small_primes_lookup.size(); ++i)
				cleared += (1.0 - cleared) * (1.0 / small_primes_lookup[i]);
			return size_t((1.0 - cleared) * sieve.size() * sieve_steps * 2);
		}();
		static size_t candidates[candidates_capacity]{}; // Intellisense false-positive

		// Start the clock after setup
		const auto start = current_time_in_ms();



		count_passes(std::cout << "(counting passes)\n");
		count_passes(size_t a, b, c, d, e, f, g, h, i, j, k, l, m, passes);
		count_passes(a = b = c = d = e = f = g = h = i = j = k = l = m = passes = 0);

		// (condition optimizes out when not benchmarking)
		while (benchmark_mode ? number < bm_stop : true)
		{
			const size_t number_before_tests = number;

			size_t* candidates_end = candidates;

			for (size_t sieve_step = 0; sieve_step < sieve_steps; ++sieve_step)
			{
				// Perform additional sieving on the static sieve
				util::vectorized_copy_n<sieve.size()>((uint256_t*)sieve.data(), (uint256_t*)static_sieve.data());
				partial_sieve(sieve);

				// Collect candidates that have not been marked composite by the sieve
				candidates_end = gather_sieve_results(candidates_end, sieve.data(), sieve.data() + sieve.size(),
													  number + (sieve_step * sieve.size() * 2));
			}
			count_passes(a += (candidates_end - candidates));



			// Collect candidates that have a prime number of bits set
			candidates_end = prime_popcount_test(candidates, candidates_end);
			count_passes(b += (candidates_end - candidates));

			// Collect candidates with an alternating bitsum that shares a GCD of 1 with a product of primes
			candidates_end = gcd_test(candidates, candidates_end);
			count_passes(c += (candidates_end - candidates));



			// Perform some div tests separately to remove some of the branchiest branches

			// base 3 mod 5, base 4 mod 17, and bases 5 and 8 mod 13 (4 remainders)
			candidates_end = div_tests_with_four_rems(candidates, candidates_end);
			count_passes(d += (candidates_end - candidates));

			// base 4 mod 7, and bases 3 and 9 mod 13 (3 remainders)
			candidates_end = div_tests_with_three_rems(candidates, candidates_end);
			count_passes(e += (candidates_end - candidates));

			// bases 3 and 5 mod 7, and 4 and 10 mod 13 (6 remainders)
			candidates_end = div_tests_with_six_rems(candidates, candidates_end);
			count_passes(f += (candidates_end - candidates));

			// bases 3, 4, 5, and 9 mod 11 (5 remainders)
			candidates_end = div_tests_with_five_rems(candidates, candidates_end);
			count_passes(g += (candidates_end - candidates));

			// bases 6, 7, and 8 mod 11 (10 remainders)
			candidates_end = div_tests_with_10_rems(candidates, candidates_end);
			count_passes(h += (candidates_end - candidates));

			// bases 6, 7, and 11 mod 13 (12 remainders)
			candidates_end = div_tests_with_12_rems(candidates, candidates_end);
			count_passes(i += (candidates_end - candidates));

			// bases 3, 5, 6, 7, 10, 11 and 12 mod 17 (16 remainders)
			candidates_end = div_tests_with_16_rems(candidates, candidates_end);
			count_passes(j += (candidates_end - candidates));



			// Check for small prime factors across all bases
			candidates_end = branchless_div_tests(candidates, candidates_end, div_test::div_tests.data(), 5);
			count_passes(k += (candidates_end - candidates));
			candidates_end = multibase_div_tests(candidates, candidates_end, div_test::div_tests.data() + 5);
			count_passes(l += (candidates_end - candidates));



			// Do full primality tests, starting with base 2
			for (size_t* candidate = candidates; candidate < candidates_end; )
			{
				number = *candidate++;

				if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

				count_passes(++m);

				// convert uint64_t to char array of ['0', '1'...] for MPIR
				char bin_str[64 + 1];
				auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
				*result.ptr = '\0';

				mpz_number.set_str(bin_str, 3);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 4);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 5);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 6);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 7);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 8);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 9);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 8); continue; }

				mpz_number.set_str(bin_str, 10);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 9); continue; }

				mpz_number.set_str(bin_str, 11);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 10); continue; }

				mpz_number.set_str(bin_str, 12);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 11); continue; }

				log_result(number, 12);
			}

			number = number_before_tests + (2 * sieve.size() * sieve_steps);



			if (next_div_test_reorder <= number)
			{
				update_div_test_order();
				next_div_test_reorder += div_test::reorder_interval;

			#if analyze_div_tests
				static auto last_perf_time = current_time_in_ms();
				static size_t last_n = 0;

				std::stringstream ss;
				ss << (number - bm_start) / 1'000'000'000 << " B ints searched.";
				const auto perf_time = current_time_in_ms();
				if (last_n != 0)
				{
					const double elapsed_seconds = double(perf_time - last_perf_time) / 1'000.0;
					const double ints_per_second = double(number - last_n) / elapsed_seconds;
					ss << ' ' << size_t((ints_per_second / 1'000'000.0) + 0.5) << " M ints/second";
				}
				ss << '\n';
				std::cout << ss.str();

				last_perf_time = perf_time;
				last_n = number;
			#endif
			}

			count_passes(++passes);
		} // end main loop



	#if analyze_div_tests
		// Run one final step before exiting
		// run_div_test_analysis();
		print_div_tests();
	#endif

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";

		log_pass_counts("Passed sieve:          ", a, (bm_size / 2));
		log_pass_counts("Passed popcount test:  ", b, a);
		log_pass_counts("Passed GCD test:       ", c, b);
		log_pass_counts("Passed 4-rem tests:    ", d, c);
		log_pass_counts("Passed 3-rem tests:    ", e, d);
		log_pass_counts("Passed 6-rem tests:    ", f, e);
		log_pass_counts("Passed 5-rem tests:    ", g, f);
		log_pass_counts("Passed 10-rem tests:   ", h, g);
		log_pass_counts("Passed 12-rem tests:   ", i, h);
		log_pass_counts("Passed 16-rem tests:   ", j, i);
		log_pass_counts("P. branchless divtests:", k, j);
		log_pass_counts("P. branching divtests: ", l, k);
		log_pass_counts("Passed b2 BPSW test:   ", m, l);
	}

} // namespace mbp



int main()
{
	mbp::print_preamble();

	mbp::find_multibase_primes();
}
