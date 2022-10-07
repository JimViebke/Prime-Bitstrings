
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

#define VCL_NAMESPACE vcl
#define MAX_VECTOR_SIZE 256
#include "../lib/vcl/vectorclass.h"

#include "bit_pattern_tests.hpp"
#include "config.hpp"
#include "full_div_tests.hpp"
#include "hardcoded_div_tests.hpp"
#include "io/io.hpp"
#include "math/franken_mpir.hpp"
#include "math/math.hpp"
#include "sieve.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/sandbox.hpp"
#include "util/simd.hpp"
#include "util/types.hpp"
#include "util/utility.hpp"


namespace mbp
{
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

		static std::stringstream moves_history;
		moves_history << ' ' << moved;

		ss << moved << " div tests changed position\n";
		ss << '(' << moves_history.str() << ")\n\n";



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



	std::vector<uint8_t> build_popcounts_lookup()
	{
		std::vector<uint8_t> popcounts;
		popcounts.reserve(1ull << 16);

		for (size_t i = 0; i < (1ull << 16); ++i)
		{
			popcounts.push_back(uint8_t(pop_count(i)));
		}

		return popcounts;
	}
	const std::vector<uint8_t> pc_lookup = build_popcounts_lookup();



	static alignas(64) sieve_container sieve;

	// 2x the expected number of candidates from the sieve passes
	constexpr size_t candidates_capacity = [] {
		double cleared = 0.0;
		for (size_t i = 1; i < small_primes_lookup.size(); ++i)
			cleared += (1.0 - cleared) * (1.0 / small_primes_lookup[i]);
		return 2 * size_t((1.0 - cleared) * sieve.size() * sieve_steps);
	}();
	static alignas(64) std::array<size_t, candidates_capacity> candidates_storage;



	class find_multibase_primes
	{
	public:
		find_multibase_primes()
		{
			gmp_rand.seed(mpir_ui{ 0xdeadbeef });

			count_passes(std::cout << "(counting passes)\n");
			count_passes(a = b = c = d = e = f = g = h = i = j = k = l = m = passes = 0);
		}

		void run()
		{
			constexpr size_t loop_size = 2ull * sieve.size() * sieve_steps;

			size_t number = benchmark_mode ? bm_start : load_from_results();

			// Round starting number down to the nearest odd multiple of the sieve size
			number -= sieve.size(); // n -= k
			number -= number % (2 * sieve.size()); // n -= n % 2k
			number += sieve.size(); // n += k

			set_up_sieve_offsets_cache(number);

			size_t next_div_test_reorder = number + div_test::reorder_interval;

			// Start the clock after setup
			const auto start = util::current_time_in_ms();

			// (condition should optimize out)
			while (benchmark_mode ? number < bm_stop : true)
			{
				// The upper 32 bits of a 64 bit integer only change every 4 billion ints.
				// Detect iterations where the upper bits can not change, and allow
				// functions to optimize based on this.
				if (util::upper_32_bits_match(number, number + loop_size))
				{
					main_loop<true>(number);
				}
				else
				{
					main_loop<false>(number);
				}

				number += loop_size;



				if (next_div_test_reorder <= number)
				{
					update_div_test_order();
					next_div_test_reorder += div_test::reorder_interval;

				#if analyze_div_tests
					run_div_test_analysis();
				#endif
				}

				count_passes(++passes);
			} // end main loop



		#if analyze_div_tests
			// Run one final step before exiting
			// run_div_test_analysis();
			print_div_tests();
		#endif

			std::cout << "Finished. " << util::current_time_in_ms() - start << " ms elapsed\n";

			log_pass_counts("Passed static sieve and\n"\
							"  popcount filter:     ", a, (bm_size / 2));
			log_pass_counts("Passed sieve:          ", b, a);
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

	private:
		__forceinline void merge_static_sieve_with_popcounts(
			uint256_t* const out,
			const uint256_t* const static_sieve_in,
			const uint256_t* const pc_in,
			const uint256_t& valid_pcs)
		{
			uint256_t pcs0 = _mm256_load_si256(pc_in + 0);
			uint256_t pcs1 = _mm256_load_si256(pc_in + 1);
			uint256_t pcs2 = _mm256_load_si256(pc_in + 2);
			uint256_t pcs3 = _mm256_load_si256(pc_in + 3);
			uint256_t ymm0 = _mm256_load_si256(static_sieve_in + 0);
			uint256_t ymm1 = _mm256_load_si256(static_sieve_in + 1);
			uint256_t ymm2 = _mm256_load_si256(static_sieve_in + 2);
			uint256_t ymm3 = _mm256_load_si256(static_sieve_in + 3);

			pcs0 = _mm256_shuffle_epi8(valid_pcs, pcs0);
			pcs1 = _mm256_shuffle_epi8(valid_pcs, pcs1);
			pcs2 = _mm256_shuffle_epi8(valid_pcs, pcs2);
			pcs3 = _mm256_shuffle_epi8(valid_pcs, pcs3);
			ymm0 = _mm256_and_si256(ymm0, pcs0);
			ymm1 = _mm256_and_si256(ymm1, pcs1);
			ymm2 = _mm256_and_si256(ymm2, pcs2);
			ymm3 = _mm256_and_si256(ymm3, pcs3);

			_mm256_store_si256(out + 0, ymm0);
			_mm256_store_si256(out + 1, ymm1);
			_mm256_store_si256(out + 2, ymm2);
			_mm256_store_si256(out + 3, ymm3);
		}

		uint256_t manually_generate_popcounts(const size_t number)
		{
			uint256_t pcs{};
			for (size_t idx = 0; idx < 32; ++idx)
			{
				pcs.m256i_u8[idx] = (build_tiny_primes_lookup() >> pop_count(number + (2 * idx))) & 1u;
			}
			return pcs;
		}

		// Out and in must be aligned on a multiple of 32 bytes.
		// Size must be at least 4 * 32 bytes.
		template<size_t n_bytes>
		__forceinline void copy_static_sieve_with_popcount_filter(uint256_t* out,
																  const uint256_t* in,
																  size_t number)
		{
			constexpr size_t high_bits_mask = size_t(-1) << (16 + 1); // 64 bits == 47 high bits + (16 + 1)
			constexpr size_t low_17_bits_mask = ~high_bits_mask;
			constexpr size_t prime_lookup = build_tiny_primes_lookup() >> 1; // preshift by 1 because the bottom bit is always set

			constexpr size_t leftover_bytes = n_bytes % (32 * 4);
			const uint256_t* const aligned_end = in + ((n_bytes - leftover_bytes) / 32);

			// while in != aligned_end:
			// - iterate by 4*32 blocks until we hit aligned_end, or a rollover of bits 1-16, whichever happens first
			// - if we broke before aligned_end, we hit a rollover
			// - - handle the next 4*32 blocks manually, then continue
			// handle end stuff

			for (; in != aligned_end; )
			{
				// generate a 16+16 element lookup to mark which popcounts of bits 1-16 add up to a valid prime
				const size_t high_bits = number & high_bits_mask;
				const uint16_t shifted_prime_lookup = uint16_t(prime_lookup >> pop_count(high_bits));
				const uint256_t valid_pcs = util::expand_16_bits_to_bytes(shifted_prime_lookup);

				// set up the popcount lookup pointer
				const size_t low_17_bits = number & low_17_bits_mask;
				const uint256_t* pc_lookup_ptr = (uint256_t*)&pc_lookup[low_17_bits >> 1];

				/*
				We are iterating 4*32 elements at a time, therefore we are reading 4*32 elements
				from the popcount lookup per iteration. How many times can we do this before
				reaching its end?

				AVX2 shuffle can handle 16 unique values, but a popcount of 0-16 is 17 values.
				Resolve this by treating 16 (popcount of 0xFFFF) as a special case.

				Calculate this up front so the hot loop can just run until "stop".
				*/
				size_t elements_to_FFFF = 0xFFFF - (low_17_bits >> 1);
				const size_t iters_to_FFFF = elements_to_FFFF / (4ull * 32);

				const size_t iters_to_aligned_end = (aligned_end - in) / 4; // 4 reads of 32b each

				const size_t n_iters = std::min(iters_to_FFFF, iters_to_aligned_end);
				const uint256_t* const stop = in + (n_iters * 4); // *4 because ptr resolution is 32b, not 32*4

				// process blocks of 4*32 bytes
				for (; in != stop; in += 4, out += 4, pc_lookup_ptr += 4)
				{
					merge_static_sieve_with_popcounts(out, in, pc_lookup_ptr, valid_pcs);
				}

				number += n_iters * (32ull * 4 * 2);

				// if we are not at aligned_end, we are at a rollover
				// (an 0xFFFF << 1 occurs in the next 4*32 bytes)
				if (in != aligned_end)
				{
					elements_to_FFFF -= n_iters * 32 * 4;

					// elements_to_FFFF should be less than 128.
					// >= 128 would mean that the 0xFFFF case does not occur in this 4*32 block

					// handle 4 32-bit blocks
					for (size_t block = 0; block < 4; ++block, number += 32ull * 2)
					{
						uint256_t pcs{};

						if (elements_to_FFFF >= 32) // copy 32 more bytes using existing valid_pcs
						{
							elements_to_FFFF -= 32;
							pcs = _mm256_load_si256(pc_lookup_ptr + block);
							pcs = _mm256_shuffle_epi8(valid_pcs, pcs);
						}
						else // handle 32 elements manually
						{
							pcs = manually_generate_popcounts(number);
						}

						uint256_t ymm = _mm256_load_si256(in + block);
						ymm = _mm256_and_si256(ymm, pcs);
						_mm256_store_si256(out + block, ymm);
					} // end block loop

					in += 4;
					out += 4;
					// number is already incremented
					// pc_lookup_ptr will be recalculated in the next iteration

				} // end if (in != aligned_end)
			} // end main copy/merge loop

			// <32*4 bytes left. Generate instructions for 0-3 aligned copies.

			if constexpr (leftover_bytes >= 32 * 1)
			{
				const uint256_t pcs = manually_generate_popcounts(number);
				const uint256_t ymm = _mm256_load_si256(in++);
				_mm256_store_si256(out++, _mm256_and_si256(ymm, pcs));
				number += 32ull * 2;
			}
			if constexpr (leftover_bytes >= 32 * 2)
			{
				const uint256_t pcs = manually_generate_popcounts(number);
				const uint256_t ymm = _mm256_load_si256(in++);
				_mm256_store_si256(out++, _mm256_and_si256(ymm, pcs));
				number += 32ull * 2;
			}
			if constexpr (leftover_bytes >= 32 * 3)
			{
				const uint256_t pcs = manually_generate_popcounts(number);
				const uint256_t ymm = _mm256_load_si256(in++);
				_mm256_store_si256(out++, _mm256_and_si256(ymm, pcs));
				number += 32ull * 2;
			}

			// <32 bytes left. Generate instructions for 0-1 unaligned copies.

			constexpr size_t adjust = 32 - (leftover_bytes % 32);
			if constexpr (adjust > 0)
			{
				in = (uint256_t*)(((uint8_t*)in) - adjust);
				out = (uint256_t*)(((uint8_t*)out) - adjust);

				const uint256_t pcs = manually_generate_popcounts(number - (adjust * 2));
				const uint256_t ymm = _mm256_loadu_si256(in);
				_mm256_storeu_si256(out, _mm256_and_si256(ymm, pcs));
			}
		}

		template<bool on_fast_path>
		void main_loop(const size_t number)
		{
			size_t* const candidates = candidates_storage.data();

			size_t* candidates_end = candidates;
			for (size_t sieve_step = 0; sieve_step < sieve_steps; ++sieve_step)
			{
				// Merge the static sieve and popcount data
				copy_static_sieve_with_popcount_filter<sieve.size()>(
					(uint256_t*)sieve.data(),
					(uint256_t*)static_sieve.data(),
					number + (sieve_step * sieve.size() * 2));
				count_passes(a += util::vector_count_ones(sieve.data(), sieve.size()));

				partial_sieve(sieve);

				// Collect candidates that have not been marked composite by the sieve
				candidates_end = gather_sieve_results(candidates_end, sieve.data(),
													  sieve.data() + sieve.size(),
													  number + (sieve_step * sieve.size() * 2));
			}
			count_passes(b += (candidates_end - candidates));



			// Collect candidates that have a prime number of bits set
			//candidates_end = prime_popcount_test(candidates, candidates_end);
			//count_passes(b += (candidates_end - candidates));

			// Collect candidates with an alternating bitsum that shares a GCD of 1 with a product of primes
			candidates_end = gcd_test(candidates, candidates_end);
			count_passes(c += (candidates_end - candidates));



			// Perform some div tests separately to remove some of the branchiest branches

			// base 3 mod 5, base 4 mod 17, and bases 5 and 8 mod 13 (4 remainders)
			candidates_end = div_tests_with_four_rems<on_fast_path>(candidates, candidates_end);
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
			candidates_end = div_tests_with_10_rems<on_fast_path>(candidates, candidates_end);
			count_passes(h += (candidates_end - candidates));

			// bases 6, 7, and 11 mod 13 (12 remainders)
			candidates_end = div_tests_with_12_rems<on_fast_path>(candidates, candidates_end);
			count_passes(i += (candidates_end - candidates));

			// bases 3, 5, 6, 7, 10, 11 and 12 mod 17 (16 remainders)
			candidates_end = div_tests_with_16_rems<on_fast_path>(candidates, candidates_end);
			count_passes(j += (candidates_end - candidates));



			// Check for small prime factors across all bases
			candidates_end = branchless_div_tests<on_fast_path>(candidates, candidates_end, div_test::div_tests.data(), 5);
			count_passes(k += (candidates_end - candidates));

			candidates_end = multibase_div_tests<on_fast_path>(candidates, candidates_end, div_test::div_tests.data() + 5);
			count_passes(l += (candidates_end - candidates));



			// Do full primality tests, starting with base 2
			for (size_t* candidate_ptr = candidates; candidate_ptr < candidates_end; ++candidate_ptr)
			{
				const size_t candidate = *candidate_ptr;

				if (!franken::mpir_is_likely_prime_BPSW(candidate)) continue;

				count_passes(++m);

				// convert uint64_t to char array of ['0', '1'...] for MPIR
				char bin_str[64 + 1];
				auto result = std::to_chars(&bin_str[0], &bin_str[64], candidate, 2);
				*result.ptr = '\0';

				mpz_number.set_str(bin_str, 3);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 4);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 5);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 6);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 7);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 8);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 9);
				if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 8); continue; }

				mpz_number.set_str(bin_str, 10);
				if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 9); continue; }

				mpz_number.set_str(bin_str, 11);
				if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 10); continue; }

				mpz_number.set_str(bin_str, 12);
				if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 11); continue; }

				log_result(candidate, 12);
			}
		}

	private:
		gmp_randclass gmp_rand{ gmp_randinit_mt };
		mpz_class mpz_number = 0ull;

		count_passes(size_t a, b, c, d, e, f, g, h, i, j, k, l, m, passes);
	};

} // namespace mbp



int main()
{
	mbp::print_preamble();

	mbp::find_multibase_primes mbp;
	mbp.run();
}
