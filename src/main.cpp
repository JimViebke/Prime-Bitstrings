
#define VCL_NAMESPACE vcl
#include "../lib/vcl/vectorclass.h"

#include <bitset>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include "config.hpp"
#include "io/io.hpp"
#include "math/franken_mpir.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/find.hpp"
#include "util/sandbox.hpp"
#include "util/types.hpp"
#include "util/utility.hpp"


namespace mbp
{
	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(),
														 static_sieve_primes.end(), size_t(1), std::multiplies());

	using sieve_t = uint8_t;
	// using sieve_container = std::vector<sieve_t>;
	using sieve_container = std::array<sieve_t, static_sieve_size>;
	const sieve_container generate_static_sieve()
	{
		// sieve_container sieve(static_sieve_size, true);
		sieve_container sieve{}; std::fill(begin(sieve), end(sieve), true);

		// for each prime, mark off all multiples
		for (const auto p : static_sieve_primes)
			for (size_t i = 0; i < sieve.size(); i += p)
				sieve[i] = false;

		return sieve;
	}
	static alignas(sieve_alignment) const sieve_container static_sieve = generate_static_sieve();

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

	tests_are_inlined void partial_sieve(sieve_container& sieve)
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
				// skip +0
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
		for (; sieve_ptr < sieve_end; sieve_ptr += 15, number += 30)
		{
			// By working in iterations of 15, we know every 3rd and every 5th value is false.
			// Don't check those ones.
			// 
			// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   <-- offset
			//    x  x     x        x  x        x     x  x   <-- values to check
			// x        x     x  x        x  x     x         <-- values to ignore

			*sieve_candidates = number + 2; // number + offset*2
			sieve_candidates += *(sieve_ptr + 1);

			*sieve_candidates = number + 4;
			sieve_candidates += *(sieve_ptr + 2);

			*sieve_candidates = number + 8;
			sieve_candidates += *(sieve_ptr + 4);

			*sieve_candidates = number + 14;
			sieve_candidates += *(sieve_ptr + 7);

			*sieve_candidates = number + 16;
			sieve_candidates += *(sieve_ptr + 8);

			*sieve_candidates = number + 22;
			sieve_candidates += *(sieve_ptr + 11);

			*sieve_candidates = number + 26;
			sieve_candidates += *(sieve_ptr + 13);

			*sieve_candidates = number + 28;
			sieve_candidates += *(sieve_ptr + 14);
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
			const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
			*passed_gcd_test_ptr = n0;
			if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++passed_gcd_test_ptr;

			const size_t n1 = *(candidate_ptr + 1);
			const size_t n1_pca = pop_count(n1 & 0xAAAAAAAAAAAAAAAA);
			const size_t n1_pcb = pop_count(n1 & 0x5555555555555555);
			*passed_gcd_test_ptr = n1;
			if (gcd_lookup & (1ull << (n1_pca + 32 - n1_pcb))) ++passed_gcd_test_ptr;

			const size_t n2 = *(candidate_ptr + 2);
			const size_t n2_pca = pop_count(n2 & 0xAAAAAAAAAAAAAAAA);
			const size_t n2_pcb = pop_count(n2 & 0x5555555555555555);
			*passed_gcd_test_ptr = n2;
			if (gcd_lookup & (1ull << (n2_pca + 32 - n2_pcb))) ++passed_gcd_test_ptr;

			const size_t n3 = *(candidate_ptr + 3);
			const size_t n3_pca = pop_count(n3 & 0xAAAAAAAAAAAAAAAA);
			const size_t n3_pcb = pop_count(n3 & 0x5555555555555555);
			*passed_gcd_test_ptr = n3;
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

	consteval auto generate_bitmask_lookup()
	{
		// +1 so bitmasks[n_of_rems] is always safe
		std::array<size_t, div_test::max_remainders + 1> bitmasks = { 0 };

		for (size_t i = 1; i < bitmasks.size(); ++i)
		{
			size_t bitmask = 0;
			for (size_t j = 0; j < 64; j += i)
			{
				bitmask <<= i;
				bitmask |= 1;
			}
			bitmasks[i] = bitmask;
		}

		return bitmasks;
	}
	constexpr std::array<size_t, div_test::max_remainders + 1> bitmask_lookup = generate_bitmask_lookup();

	// takes N^2 memory, even though we only need (N^2) / 2
	using popcount_t = uint16_t;
	static std::array<mbp::aligned64<popcount_t, 64>, 64> popcounts{};

	tests_are_inlined bool has_small_divisor(const size_t number)
	{
		using namespace div_test;

		if (recursive_is_divisible_by<5, in_base<3>>(number)) return true;

		if (recursive_is_divisible_by<7, in_base<3>>(number)) return true;
		if (recursive_is_divisible_by<7, in_base<4>>(number)) return true;
		if (recursive_is_divisible_by<7, in_base<5>>(number)) return true;

	#if analyze_div_tests
		bool found_div = false;
	#endif

		bool which_way_boss = div_tests.front().is_first_with_n_remainders;
		n_of_remainders_t n_of_rems_boss = div_tests.front().n_of_remainders;

		for (div_test_const auto& div_test : div_tests)
		{
			size_t rem = 0;

			//const size_t n_of_rems = div_test.n_of_remainders;
			//__assume(n_of_rems > 0);
			//__assume(n_of_rems <= max_remainders);

			const auto& my_rems = div_test.remainders;

			__assume(n_of_rems_boss > 0);
			__assume(n_of_rems_boss <= max_remainders);


			if (which_way_boss)
			{
				which_way_boss = *((&div_test.is_first_with_n_remainders) + sizeof(div_test));

				const size_t my_bitmask = bitmask_lookup[n_of_rems_boss];
				auto& my_pcs = popcounts[n_of_rems_boss];

				// for switch (n), run cases n through 1, where the index is n-1 through 0
				constexpr size_t start = __LINE__ + 10;
			#define IDX(n) ((max_remainders - (n - start)) - 1)
			#define CASE(n) [[fallthrough]]; case(IDX(n) + 1): \
				{ \
				const auto pc = pop_count(number & (my_bitmask << IDX(n))); \
					my_pcs[IDX(n)] = popcount_t(pc); \
					rem += pc * my_rems[IDX(n)]; \
				}
				switch (n_of_rems_boss) // handle cases N through 1
				{
					CASE(__LINE__); // case (max)
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__); // case (1)
					static_assert(start + max_remainders == __LINE__);
					break;
				default:
					__assume(false);
				}
			#undef CASE
			#undef IDX
			}
			else
			{
				which_way_boss = *((&div_test.is_first_with_n_remainders) + sizeof(div_test));

				const auto& my_pcs = popcounts[n_of_rems_boss];

				// for switch (n), run cases n through 1, where the index is n-1 through 0
				constexpr size_t start = __LINE__ + 5;
			#define IDX(n) ((max_remainders - (n - start)) - 1)
			#define CASE(n) [[fallthrough]]; case(IDX(n) + 1): rem += size_t(my_pcs[IDX(n)]) * my_rems[IDX(n)];
				switch (n_of_rems_boss)
				{
					CASE(__LINE__); // case (max)
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__); // case (1)
					static_assert(start + max_remainders == __LINE__);
					break;
				default:
					__assume(false);
				}
			#undef CASE
			#undef IDX
			}

			n_of_rems_boss = *((&div_test.n_of_remainders) + sizeof(div_test));

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

		auto w = std::setw;
		for (const auto& dt : div_tests)
		{
			std::cout << "   base " << std::setfill(' ') << w(2) << size_t(dt.base) << " % " << w(3) << size_t(small_primes_lookup[dt.prime_idx]) << ":  ";
			if (dt.hits == 0)
			{
				std::cout << "       -       ";
			}
			else
			{
				std::cout << w(8) << dt.hits << " hits  ";
			}

			std::cout << w(2) << size_t(dt.n_of_remainders) << " remainders: 1";
			for (size_t j = 1; j < dt.n_of_remainders; ++j)
			{
				std::cout << ' ' << w(3) << dt.remainders[j];
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

			// clear "firsts"
			for (auto& dt : div_tests)
			{
				dt.is_first_with_n_remainders = false;
			}

			// set firsts
			for (size_t i = 0; i <= div_test::max_remainders; ++i)
			{
				for (auto& div_test : div_tests)
				{
					if (div_test.n_of_remainders == i)
					{
						div_test.is_first_with_n_remainders = true;
						break;
					}
				}
			}

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

		static alignas(sieve_alignment) sieve_container sieve = static_sieve;

		// Round starting number down to the nearest odd multiple of the sieve sieze
		number -= static_sieve.size(); // n -= k
		number -= number % (2 * static_sieve.size()); // n -= n % 2k
		number += static_sieve.size(); // n += k

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
			return size_t((1.0 - cleared) * double(static_sieve_size) * 1.1);
		}();
		static size_t scratch[scratch_size]{};

		// Start the clock after setup
		const auto start = current_time_in_ms();



		count_passes(size_t a, b, c, d);
		count_passes(a = b = c = d = 0);

		// (condition optimizes out when not benchmarking)
		while (benchmark_mode ? number < bm_stop : true)
		{
			// Perform additional sieving on the static sieve
			vectorized_copy((__m256i*) sieve.data(),
							(__m256i*) static_sieve.data(),
							static_sieve.size());
			partial_sieve(sieve);

			const size_t number_before_tests = number;



			// 1. Collect candidates that have not been marked composite by the sieve
			const size_t* const sieve_candidates = gather_sieve_results(scratch, sieve.data(), sieve.data() + sieve.size(), number);

			count_passes(a += (sieve_candidates - scratch));

			// 2. Collect candidates that have a prime number of bits set
			const size_t* const passed_pc_test_ptr = prime_popcount_test(scratch, sieve_candidates, tiny_primes_lookup);

			count_passes(b += (passed_pc_test_ptr - scratch));

			// 3. Collect candidates with an alternating bitsum that shares a GCD of 1 with a product of primes
			const size_t* const passed_gcd_test_ptr = gcd_test(scratch, passed_pc_test_ptr, gcd_lookup);

			count_passes(c += (passed_gcd_test_ptr - scratch));

			for (size_t* n = scratch; n < passed_gcd_test_ptr; ++n)
			{
				number = *n;

				// Bail if n has a small prime factor in any base
				if (has_small_divisor(number)) continue;

				count_passes(++d);



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
			for (const auto& dt : div_test::div_tests)
			{
				if (dt.hits >= next_div_test_checkpoint)
				{
					run_div_test_analysis();
					next_div_test_checkpoint += div_test_log_interval;
					break;
				}
			}
		#endif

		} // end main loop

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";

		count_passes(std::cout << "Passed sieve:      " << std::setw(9) << a << '\n');
		count_passes(std::cout << "Passed prime test: " << std::setw(9) << b << '\n');
		count_passes(std::cout << "Passed GCD test:   " << std::setw(9) << c << '\n');
		count_passes(std::cout << "Passed div tests:  " << std::setw(9) << d << '\n');
	}

} // namespace mbp



int main()
{
	mbp::print_preamble();

	mbp::find_multibase_primes();
}
