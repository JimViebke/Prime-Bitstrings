
// Hacky way to hide the contents of zmmintrin.h
// For some reason MSVC always includes it from immintrin.h
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
	// compile time, based on requiring primes to make at least (sieve_size / X) number of writes, where X
	// is a tuneable knob.
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
		sieve_prime_t* cache_ptr = sieve_offsets_cache.data() + static_sieve_primes.size() + 1;

		for (;;)
		{
			// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
			//    x  x     x        x  x        x     x  x   <-- values to mark composite/false
			// x        x     x  x        x  x     x         <-- values to ignore

			// Get the next prime
			const size_t p = *prime_ptr;

			const size_t fifteen_p = p * 15;

			// Get the index of the next odd multiple of p*3
			sieve_t* j = sieve_begin + *cache_ptr;

			*(((j - (14 * p)) >= sieve_begin) ? (j - (14 * p)) : j) = false;
			*(((j - (13 * p)) >= sieve_begin) ? (j - (13 * p)) : j) = false;
			*(((j - (11 * p)) >= sieve_begin) ? (j - (11 * p)) : j) = false;
			*(((j - (8 * p)) >= sieve_begin) ? (j - (8 * p)) : j) = false;
			*(((j - (7 * p)) >= sieve_begin) ? (j - (7 * p)) : j) = false;
			*(((j - (4 * p)) >= sieve_begin) ? (j - (4 * p)) : j) = false;
			*(((j - (2 * p)) >= sieve_begin) ? (j - (2 * p)) : j) = false;
			*(((j - p) >= sieve_begin) ? (j - p) : j) = false;

			// Stop marking 15p early, then handle cleanup after the loop.
			const sieve_t* const padded_end = sieve_end - fifteen_p;

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

				j += fifteen_p;
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
			*cache_ptr = sieve_prime_t((j + fifteen_p) - sieve_end);

			++prime_ptr;
			++cache_ptr;

			if (p == last_prime_for_stepping_by_fifteen) break;
		}

		for (;;)
		{
			// Get the next prime
			const size_t p = *prime_ptr;
			// For readability:
			const size_t two_p = p * 2;
			const size_t three_p = p * 3;

			// Get the index of the next odd multiple of p*3
			sieve_t* j = sieve_begin + *cache_ptr;

			// j is aligned to a multiple of 3, so we have 0, 1, or 2 multiples of p *behind* j that must be crossed off.
			// Use two branchless writes to j-p and j-2p if they exist, falling back to a write to j in both cases
			*(((j - p) >= sieve_begin) ? (j - p) : j) = false;
			*(((j - two_p) >= sieve_begin) ? (j - two_p) : j) = false;

			// Stop marking 3p early, then handle cleanup after the loop.
			const sieve_t* const padded_end = sieve_end - three_p;

			// Each iteration, j points to an (already marked) multiple of p and 3.
			// Only mark j + p and j + 2p.
			do
			{
				*(j + p) = false;
				*(j + two_p) = false;
				j += three_p;
			} while (j < padded_end);

			// Same as above, we have 0, 1, or 2 remaining writes to perform
			*(((j + p) < sieve_end) ? (j + p) : j) = false;
			*(((j + two_p) < sieve_end) ? (j + two_p) : j) = false;

			// Calculate the offset and update the cache for the next sieving
			*cache_ptr = sieve_prime_t((j + three_p) - sieve_end);

			++prime_ptr;
			++cache_ptr;

			if (p == last_prime_for_stepping_by_threes) break;
		}

		// Sieve the remaining primes, which are likely too large relative to the sieve size
		// to benefit from stepping by more than p.

		for (; prime_ptr < small_primes_lookup.data() + small_primes_lookup.size();
			 ++prime_ptr, ++cache_ptr)
		{
			// Get the next prime
			const size_t p = *prime_ptr;

			// Get the index of the next odd multiple of p
			sieve_t* j = sieve_begin + *cache_ptr;

			// Mark off each (implicitly) odd multiple of p
			do
			{
				*j = false;
				j += p;
			} while (j < sieve_end);

			// Calculate the offset and update the cache for the next sieving
			*cache_ptr = sieve_prime_t(j - sieve_end);
		}
	}



	tests_are_inlined const size_t* const gather_sieve_results(size_t* sieve_candidates,
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

		// The offset must be multiplied by 2 when calculating the bitstring, because the sieve only represents odd values.
		// There is a performance gotcha here: when this offset reaches 128, almost half of the encoded instructions below
		// become slightly larger, and no longer fit in 16 bytes. We can fix this by advancing both "number" and
		// "sieve_ptr" halfway through the loop, and then continuing with new, smaller offsets. We make these offsets smaller
		// by subtracting offset[23], the midpoint, from each of the following 24 offsets.

		constexpr std::array<size_t, 48> offsets = []() consteval {
			std::array<size_t, 48> offsets{};
			size_t idx = 0;
			// Find the offsets that are not divisible by 3, 5, or 7
			for (size_t n = 1; n <= (3ull * 5 * 7); ++n)
				if (n % 3 != 0 && n % 5 != 0 && n % 7 != 0)
					offsets[idx++] = n;

			// Perform scaling on the second half of the offsets, per above
			for (size_t n = 24; n < offsets.size(); ++n)
				offsets[n] -= offsets[23];

			return offsets;
		}();

		for (; sieve_ptr < sieve_end;
			 sieve_ptr += (3ull * 5 * 7) - offsets[23],
			 number += (3ull * 5 * 7 * 2) - (offsets[23] * 2))
		{
			auto it = offsets.cbegin();

			*sieve_candidates = number + *it * 2; // unconditionally store number + offset*2
			sieve_candidates += *(sieve_ptr + *it++); // conditionally increment
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // write #5

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 10

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 15

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 20

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);

			number += *it * 2;
			sieve_ptr += *it++;
			*sieve_candidates = number;
			sieve_candidates += *(sieve_ptr); // 24 is a special case - see above

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 25

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 30

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 35

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 40

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 45

			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++);
			*sieve_candidates = number + *it * 2;
			sieve_candidates += *(sieve_ptr + *it++); // 48
		}

		return sieve_candidates;
	}

	tests_are_inlined const size_t* const prime_popcount_test(size_t* passed_pc_test_ptr,
															  const size_t* const candidates_end,
															  const size_t& tiny_primes_lookup)
	{
		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - passed_pc_test_ptr) & 0b11);
		const size_t* candidate_ptr = passed_pc_test_ptr;

		for (; candidate_ptr < candidates_end_rounded; candidate_ptr += 4)
		{
			const size_t c0 = *candidate_ptr;
			const size_t pc_c0 = pop_count(c0);
			*passed_pc_test_ptr = c0;
			if (tiny_primes_lookup & (1ull << pc_c0)) ++passed_pc_test_ptr;

			const size_t c1 = *(candidate_ptr + 1);
			const size_t pc_c1 = pop_count(c1);
			*passed_pc_test_ptr = c1;
			if (tiny_primes_lookup & (1ull << pc_c1)) ++passed_pc_test_ptr;

			const size_t c2 = *(candidate_ptr + 2);
			const size_t pc_c2 = pop_count(c2);
			*passed_pc_test_ptr = c2;
			if (tiny_primes_lookup & (1ull << pc_c2)) ++passed_pc_test_ptr;

			const size_t c3 = *(candidate_ptr + 3);
			const size_t pc_c3 = pop_count(c3);
			*passed_pc_test_ptr = c3;
			if (tiny_primes_lookup & (1ull << pc_c3)) ++passed_pc_test_ptr;
		}

		// handle last few elements
		for (; candidate_ptr < candidates_end; ++candidate_ptr)
		{
			const size_t c0 = *candidate_ptr;
			const size_t pc_c0 = pop_count(c0);
			*passed_pc_test_ptr = c0;
			if (tiny_primes_lookup & (1ull << pc_c0)) ++passed_pc_test_ptr;
		}

		return passed_pc_test_ptr;
	}

	tests_are_inlined const size_t* const gcd_test(size_t* passed_gcd_test_ptr,
												   const size_t* const candidates_end,
												   const size_t& gcd_lookup)
	{
		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - passed_gcd_test_ptr) & 0b11);
		const size_t* candidate_ptr = passed_gcd_test_ptr;

		for (; candidate_ptr < candidates_end_rounded; candidate_ptr += 4)
		{
			const size_t n0 = *candidate_ptr;
			*passed_gcd_test_ptr = n0;
			const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++passed_gcd_test_ptr;

			const size_t n1 = *(candidate_ptr + 1);
			*passed_gcd_test_ptr = n1;
			const size_t n1_pca = pop_count(n1 & 0xAAAAAAAAAAAAAAAA);
			const size_t n1_pcb = pop_count(n1 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n1_pca + 32 - n1_pcb))) ++passed_gcd_test_ptr;

			const size_t n2 = *(candidate_ptr + 2);
			*passed_gcd_test_ptr = n2;
			const size_t n2_pca = pop_count(n2 & 0xAAAAAAAAAAAAAAAA);
			const size_t n2_pcb = pop_count(n2 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n2_pca + 32 - n2_pcb))) ++passed_gcd_test_ptr;

			const size_t n3 = *(candidate_ptr + 3);
			*passed_gcd_test_ptr = n3;
			const size_t n3_pca = pop_count(n3 & 0xAAAAAAAAAAAAAAAA);
			const size_t n3_pcb = pop_count(n3 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n3_pca + 32 - n3_pcb))) ++passed_gcd_test_ptr;
		}

		// handle last few elements
		for (; candidate_ptr < candidates_end; ++candidate_ptr)
		{
			const size_t n0 = *candidate_ptr;
			const auto n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const auto n0_pcb = pop_count(n0 & 0x5555555555555555);

			*passed_gcd_test_ptr = n0;
			passed_gcd_test_ptr += bool(gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb)));
		}

		return passed_gcd_test_ptr;
	}

	tests_are_inlined const size_t* const div_by_5_test(size_t* input,
														const size_t* const candidates_end)
	{
		using namespace div_test;

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			// Only advance the pointer if the number is still a candidate
			if (!recursive_is_divisible_by<5, in_base<3>>(number)) ++output;
		}

		return output;
	}

	tests_are_inlined const size_t* const div_by_7_base_4_test(size_t* input,
															   const size_t* const candidates_end)
	{
		using namespace div_test;

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			// Only advance the pointer if the number is still a candidate
			if (!recursive_is_divisible_by<7, in_base<4>>(number)) ++output;
		}

		return output;
	}

	tests_are_inlined const size_t* const div_by_7_bases_3_and_5_test(size_t* input,
																	  const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<3, 7>::val;
		static_assert(bitmask == bitmask_for<5, 7>::val);
		static_assert(period_of<bitmask>::val == 6);

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
			size_t b3_rem = pc_0;
			size_t b5_rem = pc_0;

			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b3_rem += pc_1 * pow_mod<3, 1, 7>::rem;
			b5_rem += pc_1 * pow_mod<5, 1, 7>::rem;

			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b3_rem += pc_2 * pow_mod<3, 2, 7>::rem;
			b5_rem += pc_2 * pow_mod<5, 2, 7>::rem;

			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b3_rem += pc_3 * pow_mod<3, 3, 7>::rem;
			b5_rem += pc_3 * pow_mod<5, 3, 7>::rem;

			const size_t pc_4 = pop_count(number & (bitmask << 4));
			b3_rem += pc_4 * pow_mod<3, 4, 7>::rem;
			b5_rem += pc_4 * pow_mod<5, 4, 7>::rem;

			const size_t pc_5 = pop_count(number & (bitmask << 5));
			b3_rem += pc_5 * pow_mod<3, 5, 7>::rem;
			b5_rem += pc_5 * pow_mod<5, 5, 7>::rem;

			// Only advance the pointer if the number is still a candidate
			bool still_a_candidate = true;
			if (has_small_prime_factor(b3_rem, get_prime_index<7>::idx)) still_a_candidate = false;
			if (has_small_prime_factor(b5_rem, get_prime_index<7>::idx)) still_a_candidate = false;
			output += still_a_candidate;
		}

		return output;
	}

	tests_are_inlined const size_t* const div_tests_b8m13_b5m13_b4m17(size_t* input,
																	  const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// base  8 % 13:   4 remainders: 1   8  12   5
		// base  5 % 13:   4 remainders: 1   5  12   8
		// base  4 % 17:   4 remainders: 1   4  16  13

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<8, 13>::val;
		static_assert(bitmask == bitmask_for<5, 13>::val &&
					  bitmask == bitmask_for<4, 17>::val);
		static_assert(period_of<bitmask>::val == 4);

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
			size_t b8m13_rem = pc_0;
			size_t b5m13_rem = pc_0;
			size_t b4m17_rem = pc_0;

			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b8m13_rem += pc_1 * pow_mod<8, 1, 13>::rem;
			b5m13_rem += pc_1 * pow_mod<5, 1, 13>::rem;
			b4m17_rem += pc_1 * pow_mod<4, 1, 17>::rem;

			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b8m13_rem += pc_2 * pow_mod<8, 2, 13>::rem;
			b5m13_rem += pc_2 * pow_mod<5, 2, 13>::rem;
			b4m17_rem += pc_2 * pow_mod<4, 2, 17>::rem;

			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b8m13_rem += pc_3 * pow_mod<8, 3, 13>::rem;
			b5m13_rem += pc_3 * pow_mod<5, 3, 13>::rem;
			b4m17_rem += pc_3 * pow_mod<4, 3, 17>::rem;

			// Only advance the pointer if the number is still a candidate
			bool still_a_candidate = true;
			if (has_small_prime_factor(b8m13_rem, get_prime_index<13>::idx)) still_a_candidate = false;
			if (has_small_prime_factor(b5m13_rem, get_prime_index<13>::idx)) still_a_candidate = false;
			if (has_small_prime_factor(b4m17_rem, get_prime_index<17>::idx)) still_a_candidate = false;
			output += still_a_candidate;
		}

		return output;
	}



	tests_are_inlined bool has_small_divisor(const size_t number)
	{
		using namespace div_test;

		//if (recursive_is_divisible_by<7, in_base<3>>(number)) return true;
		//if (recursive_is_divisible_by<7, in_base<4>>(number)) return true;
		//if (recursive_is_divisible_by<7, in_base<5>>(number)) return true;

	#if analyze_div_tests
		bool found_div = false;
	#endif

		static_assert(sizeof(remainder_t) == 1);

		// Convert each half of number to a 32-byte bitmask
		const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));
		const uint256_t mask_upper = util::expand_bits_to_bytes(number >> 32);

		uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)&div_tests[0].remainders[0]);
		uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)&div_tests[0].remainders[32]);

		for (div_test_const div_test_t& div_test : div_tests)
		{
			const uint256_t rems_lower = _mm256_and_si256(mask_lower, ymm0);
			const uint256_t rems_upper = _mm256_and_si256(mask_upper, ymm1);

			ymm0 = _mm256_loadu_si256((uint256_t*)(((uint8_t*)&div_test.remainders) + sizeof(div_test_t) + 0));
			ymm1 = _mm256_loadu_si256((uint256_t*)(((uint8_t*)&div_test.remainders) + sizeof(div_test_t) + 32));

			const size_t rem = util::vcl_hadd2_x(rems_upper, rems_lower);

			if (has_small_prime_factor(rem, div_test.prime_idx))
			{
			#if analyze_div_tests
				div_test.hits++;
				found_div = true;
				return true;
			#else
				return true;
			#endif
			}
		}

	#if analyze_div_tests
		return found_div;
	#else
		return false;
	#endif
	}

	void print_div_tests()
	{
	#if analyze_div_tests
		using namespace div_test;

		//std::sort(div_tests.begin(), div_tests.end(), [] (const auto& a, const auto& b)
		//		  {
		//			  //return a.hits < b.hits;
		//			  //return a.n_of_remainders < b.n_of_remainders;

		//			  if (a.n_of_remainders == b.n_of_remainders)
		//				  return a.base < b.base;
		//			  else
		//				  return a.n_of_remainders < b.n_of_remainders;

		//			  //if (a.prime_idx == b.prime_idx)
		//				 // return a.base < b.base;
		//			  //else
		//				 // return a.prime_idx < b.prime_idx;

		//			  //return
		//				 // a.hits / a.n_of_remainders <
		//				 // b.hits / b.n_of_remainders;
		//		  });

		std::cout << '\n';

		auto w = std::setw;
		for (const auto& dt : div_tests)
		{
			std::cout << "   base " << std::setfill(' ') << w(2) << size_t(dt.base) << " % " << w(3) << size_t(small_primes_lookup[dt.prime_idx]) << ":  ";
			//if (dt.hits == 0)
			//{
			//	std::cout << "       -       ";
			//}
			//else
			//{
			//	std::cout << w(8) << dt.hits << " hits  ";
			//}

			std::cout << w(2) << size_t(dt.n_of_remainders) << " remainders: 1";
			for (size_t j = 1; j < dt.n_of_remainders; ++j)
			{
				std::cout << ' ' << w(3) << size_t(dt.remainders[j]);
				if (j == 20)
				{
					std::cout << " ...";
					break;
				}
			}
			std::cout << '\n';
		}
	#endif
	}

	void run_div_test_analysis()
	{
	#if analyze_div_tests
		using namespace div_test;

		auto div_test_pred = [](auto a, auto b) { return a.hits < b.hits; };

		if (std::is_sorted(div_tests.rbegin(), div_tests.rend(), div_test_pred))
		{
			std::cout << "Div tests have not changed frequency ordering\n";
		}
		else
		{
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
			std::cout << '(' << ss.str() << ")\n";
		}
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

		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();
		constexpr size_t gcd_lookup = build_gcd_lookup();

	#if analyze_div_tests
		const size_t div_test_log_interval = 1'000'000;
		size_t next_div_test_checkpoint = div_test_log_interval;
	#endif

		constexpr size_t scratch_size = [] {
			double cleared = 0.0;
			for (size_t i = 1; i < small_primes_lookup.size(); ++i)
				cleared += ((1.0 - cleared) * (1.0 / small_primes_lookup[i]));
			return size_t((1.0 - cleared) * double(sieve.size()) * 1.1);
		}();
		static size_t scratch[scratch_size]{};

		// Start the clock after setup
		const auto start = current_time_in_ms();



		count_passes(size_t a, b, c, d, e, f, g, h);
		count_passes(a = b = c = d = e = f = g = h = 0);

		// (condition optimizes out when not benchmarking)
		while (benchmark_mode ? number < bm_stop : true)
		{
			const size_t number_before_tests = number;



			// Perform additional sieving on the static sieve
			util::vectorized_copy_n<sieve.size()>((uint256_t*)sieve.data(), (uint256_t*)static_sieve.data());
			partial_sieve(sieve);



			// Collect candidates that have not been marked composite by the sieve
			const size_t* candidates_end = gather_sieve_results(scratch, sieve.data(), sieve.data() + sieve.size(), number);
			count_passes(a += (candidates_end - scratch));



			// Collect candidates that have a prime number of bits set
			candidates_end = prime_popcount_test(scratch, candidates_end, tiny_primes_lookup);
			count_passes(b += (candidates_end - scratch));



			// Collect candidates with an alternating bitsum that shares a GCD of 1 with a product of primes
			candidates_end = gcd_test(scratch, candidates_end, gcd_lookup);
			count_passes(c += (candidates_end - scratch));



			// Perform some div tests separately to remove some of the branchiest branches
			candidates_end = div_by_5_test(scratch, candidates_end);
			count_passes(d += (candidates_end - scratch));

			candidates_end = div_by_7_base_4_test(scratch, candidates_end);
			count_passes(e += (candidates_end - scratch));

			candidates_end = div_by_7_bases_3_and_5_test(scratch, candidates_end);
			count_passes(f += (candidates_end - scratch));

			candidates_end = div_tests_b8m13_b5m13_b4m17(scratch, candidates_end);
			count_passes(g += (candidates_end - scratch));



			for (size_t* candidate = scratch; candidate < candidates_end; )
			{
				number = *candidate++;

				// Bail if n has a small prime factor in any base
				if (has_small_divisor(number)) continue;

				count_passes(++h);



				// Do full primality tests, starting with base 2
				if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

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

			number = number_before_tests + 2 * sieve.size();



		#if analyze_div_tests
			print_div_tests();
			std::cin.ignore();

			//for (const auto& dt : div_test::div_tests)
			//{
			//	if (dt.hits >= next_div_test_checkpoint)
			//	{
			//		run_div_test_analysis();
			//		next_div_test_checkpoint += div_test_log_interval;
			//		break;
			//	}
			//}
		#endif

		} // end main loop

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";



		count_passes(auto w = std::setw);
		count_passes(std::cout << "Passed sieve:         " << w(10) << a << " (removed ~" << w(3) << 100 - (a * 100 / bm_size) << "%)\n");
		count_passes(std::cout << "Passed prime test:    " << w(10) << b << " (removed ~" << w(3) << 100 - (b * 100 / a) << "%)\n");
		count_passes(std::cout << "Passed GCD test:      " << w(10) << c << " (removed ~" << w(3) << 100 - (c * 100 / b) << "%)\n");
		count_passes(std::cout << "Passed b4 / 5 test:   " << w(10) << d << " (removed ~" << w(3) << 100 - (d * 100 / c) << "%)\n");
		count_passes(std::cout << "Passed b4 / 7 test:   " << w(10) << e << " (removed ~" << w(3) << 100 - (e * 100 / d) << "%)\n");
		count_passes(std::cout << "Passed b3&5 / 7 test: " << w(10) << f << " (removed ~" << w(3) << 100 - (f * 100 / e) << "%)\n");
		count_passes(std::cout << "P. b8/13 b5/13 b4/17: " << w(10) << g << " (removed ~" << w(3) << 100 - (g * 100 / f) << "%)\n");
		count_passes(std::cout << "Passed div tests:     " << w(10) << h << " (removed ~" << w(3) << 100 - (h * 100 / g) << "%)\n");

	}

} // namespace mbp



int main()
{
	mbp::print_preamble();

	mbp::find_multibase_primes();
}
