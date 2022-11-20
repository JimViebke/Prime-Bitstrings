
#include <bitset>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

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

	void run_div_test_analysis(const size_t number)
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

		std::stringstream ss;
		ss << moved << " div tests changed position\n";
		ss << '(' << moves_history.str() << ")\n\n";

		static auto last_perf_time = util::current_time_in_ms();
		static size_t last_n = 0;

		ss << (number - bm_start) / 1'000'000'000 << " B ints searched.";
		const auto perf_time = util::current_time_in_ms();
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
	#else
		number; // suppress "unused formal parameter" warning
	#endif
	}



	static const sieve_container static_sieve = generate_static_sieve();
	static sieve_container sieve;

	void copy_static_sieve_with_bit_pattern_filters(size_t number)
	{
		using block_t = uint64_t;
		constexpr size_t elements_per_block = sizeof(block_t) * 8;

		// 64 bits == 47 high bits + (16 "inner" bits) + 1 low bit (always set)
		constexpr size_t bits_1_16_mask = 0xFFFFull << 1;
		constexpr size_t outer_48_bits_mask = ~bits_1_16_mask;

		constexpr size_t leftover_elements = sieve.size() % (elements_per_block * 4);

		const block_t* in = (const block_t*)static_sieve.data();
		const block_t* const aligned_end = in + ((sieve.size() - leftover_elements) / elements_per_block);

		block_t* out = (block_t*)sieve.data();

		// while in != aligned_end:
		// - iterate by 4 blocks until we hit aligned_end, or a rollover of bits 1-16, whichever happens first
		// - if we broke before aligned_end, we hit a rollover
		//   - handle the next 4 blocks seperately, then continue
		// handle cleanup

		const block_t* pc_lookup_ptr{};
		const block_t* gcd_lookup_ptr{};

		const auto set_lookup_ptrs = [&](const size_t next_number) {
			const size_t outer_bits = next_number & outer_48_bits_mask;
			const size_t bits_1_16 = (next_number & bits_1_16_mask) >> 1;

			const size_t bit_offset = bits_1_16 & 0b111;
			const size_t byte_offset = bits_1_16 / 8;

			const auto pc_idx = pop_count(outer_bits) - 2; // -2 to normalize popcount 2-48 to idx 0-46
			pc_lookup_ptr = (const block_t*)(pc_lookup[bit_offset][pc_idx].data() + byte_offset);

			const auto outer_even_pc = pop_count(outer_bits & 0x5555555555555555);
			const auto outer_odd_pc = pop_count(outer_bits & 0xAAAAAAAAAAAAAAAA);
			const auto gcd_idx = (outer_even_pc - outer_odd_pc) + 23; // +23 to normalize -23,24 to 0,47
			gcd_lookup_ptr = (const block_t*)(gcd_lookup[bit_offset][gcd_idx].data() + byte_offset);
		};

		const auto merge_one_block = [&] {
			const size_t elements_to_rollover = (pow_2_16 - ((number & bits_1_16_mask) >> 1)) % pow_2_16; // map 65,536 -> 0

			block_t mask = *in++;
			const block_t bit_patterns_mask = (*pc_lookup_ptr++) & (*gcd_lookup_ptr++);

			if (elements_to_rollover >= elements_per_block) // copy another block using lookup data
			{
				mask &= bit_patterns_mask;
			}
			else // elements_to_rollover is 0 to 63
			{
				set_lookup_ptrs(number + (elements_to_rollover * 2));
				const block_t new_pc_gcd_mask = (*pc_lookup_ptr & *gcd_lookup_ptr) << elements_to_rollover;

				const block_t select_from_old = (1ull << elements_to_rollover) - 1; // up to and including rollover
				const block_t select_from_new = ~select_from_old;

				mask &= ((bit_patterns_mask & select_from_old) |
						 (new_pc_gcd_mask & select_from_new));

				// (re)set
				set_lookup_ptrs(number + elements_per_block * 2);
			}

			*out++ = mask;

			number += elements_per_block * 2;
		};

		for (; in != aligned_end; )
		{
			set_lookup_ptrs(number);

			// We are iterating and reading 4*elements_per_block elements per iteration.
			// How many times can we do this before reaching the lookups' ends?
			// Calculate this up front so the hot loop can run using a simpler condition.

			const size_t bits_1_16 = (number & bits_1_16_mask) >> 1;
			size_t elements_to_rollover = (1ull << 16) - bits_1_16; // elements to rollover

			const size_t n_blocks = std::min(elements_to_rollover / (4ull * elements_per_block), // 4*elements_per_block elements per iteration
											 (aligned_end - in) / 4ull); // 4 reads of elements_per_block elements each

			for (size_t offset = 0; offset < n_blocks * 4; offset += 4)
			{
				const uint256_t pc_data = _mm256_loadu_si256((uint256_t*)(pc_lookup_ptr + offset));
				const uint256_t gcd_data = _mm256_loadu_si256((uint256_t*)(gcd_lookup_ptr + offset));
				const uint256_t ss_data = _mm256_loadu_si256((uint256_t*)(in + offset));

				uint256_t merged_data = _mm256_and_si256(pc_data, gcd_data);
				merged_data = _mm256_and_si256(merged_data, ss_data);

				_mm256_storeu_si256((uint256_t*)(out + offset), merged_data);
			}

			in += n_blocks * 4;
			out += n_blocks * 4;
			pc_lookup_ptr += n_blocks * 4;
			gcd_lookup_ptr += n_blocks * 4;
			number += n_blocks * 4 * elements_per_block * 2;

			// if we are not at aligned_end, we are at a rollover
			if (in != aligned_end)
			{
				// handle 4 blocks
				for (size_t block = 0; block < 4; ++block)
				{
					merge_one_block();
				}
			}
		} // end main copy/merge loop

		// Less than elements_per_block*4 elements left. Generate instructions for 0-3 copies.

		if constexpr (leftover_elements >= elements_per_block * 1)
		{
			merge_one_block();
		}
		if constexpr (leftover_elements >= elements_per_block * 2)
		{
			merge_one_block();
		}
		if constexpr (leftover_elements >= elements_per_block * 3)
		{
			merge_one_block();
		}

		constexpr size_t adjust = leftover_elements % elements_per_block;
		if constexpr (adjust > 0)
		{
			merge_one_block();
		}
	}



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
			count_passes(a = ps15 = b = c = d = e = f = g = h = i = j = 0);
			count_passes(k = l = m = n = o = p = q = passes = pc_hash = 0);
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
					print_div_tests();
					//run_div_test_analysis(number);
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
							"  bit pattern filters: ", a, (bm_size / 2));
			log_pass_counts("Passed sieve:          ", ps15, a);
			log_pass_counts("Passed 4-rem tests:    ", b, ps15);
			log_pass_counts("Passed b13m17 test:    ", c, b);
			log_pass_counts("Passed 3-rem tests:    ", d, c);
			log_pass_counts("Passed 4-rem tests, p2:", e, d);
			log_pass_counts("Passed 8-rem tests:    ", f, e);
			log_pass_counts("Passed 6-rem tests:    ", g, f);
			log_pass_counts("Passed 5-rem tests:    ", h, g);
			log_pass_counts("Passed 10-rem tests:   ", i, h);
			log_pass_counts("Passed 16-rem tests:   ", j, i);
			log_pass_counts("Passed 12-rem tests:   ", k, j);
			log_pass_counts("P. branchless divtests:", l, k);
			log_pass_counts("P. branching divtests: ", m, l);
			log_pass_counts("Passed b2 BPSW test:   ", n, m);
			log_pass_counts("Passed b3 prime test:  ", o, n);
			log_pass_counts("Passed b4 prime test:  ", p, o);
			log_pass_counts("Passed b5 prime test:  ", q, p);

			count_passes(std::cout << "\nhash of pass counts: " <<
						 std::hex << pc_hash << std::dec << '\n');
		}

	private:

		template<bool on_fast_path>
		void main_loop(const size_t number)
		{
			size_t* const candidates = candidates_storage.data();

			size_t* candidates_end = candidates;
			for (size_t sieve_step = 0; sieve_step < sieve_steps; ++sieve_step)
			{
				// Merge the static sieve, popcount, and gcd data
				copy_static_sieve_with_bit_pattern_filters(number + (sieve_step * sieve.size() * 2));
				count_passes(a += sieve.count_bits());

				partial_sieve(sieve
							  count_passes(, ps15));

				candidates_end = gather_sieve_results(
					candidates_end, sieve, number + (sieve_step * sieve.size() * 2));
			}



			// Perform some div tests separately when a specialized implementation is faster

			// base 3 mod 5, base 4 mod 17, and bases 5 and 8 mod 13 (4 remainders)
			candidates_end = div_tests_with_four_rems<on_fast_path>(candidates, candidates_end);
			count_passes(b += (candidates_end - candidates));

			// base 13 mod 17 (4 remainders, 8 candidates)
			candidates_end = base13_mod17_div_test<on_fast_path>(candidates, candidates_end);
			count_passes(c += (candidates_end - candidates));

			// base 4 mod 7, and bases 3 and 9 mod 13 (3 remainders)
			candidates_end = div_tests_with_three_rems(candidates, candidates_end);
			count_passes(d += (candidates_end - candidates));

			// base 12 mod 19, and 6 mod 37 (4 remainders, part 2)
			candidates_end = two_div_tests_with_four_rems<on_fast_path>(candidates, candidates_end);
			count_passes(e += (candidates_end - candidates));

			// bases 8 and 9 mod 17 (8 remainders)
			candidates_end = div_tests_with_8_rems<on_fast_path>(candidates, candidates_end);
			count_passes(f += (candidates_end - candidates));

			// bases 3 and 5 mod 7, and 4 and 10 mod 13 (6 remainders)
			candidates_end = div_tests_with_six_rems<on_fast_path>(candidates, candidates_end);
			count_passes(g += (candidates_end - candidates));

			// bases 3, 4, 5, and 9 mod 11 (5 remainders)
			candidates_end = div_tests_with_five_rems(candidates, candidates_end);
			count_passes(h += (candidates_end - candidates));

			// bases 6, 7, and 8 mod 11 (10 remainders)
			candidates_end = div_tests_with_10_rems<on_fast_path>(candidates, candidates_end);
			count_passes(i += (candidates_end - candidates)); // fixme

			// bases 3, 5, 6, 7, 10, 11 and 12 mod 17 (16 remainders)
			candidates_end = div_tests_with_16_rems<on_fast_path>(candidates, candidates_end);
			count_passes(j += (candidates_end - candidates)); // fixme

			// bases 6, 7, and 11 mod 13 (12 remainders)
			candidates_end = div_tests_with_12_rems<on_fast_path>(candidates, candidates_end);
			count_passes(k += (candidates_end - candidates));



			// Check for small prime factors across all bases
			candidates_end = branchless_div_tests<on_fast_path>(candidates, candidates_end, div_test::div_tests.data(), div_test::n_of_branchless_tests);
			count_passes(l += (candidates_end - candidates));

			candidates_end = multibase_div_tests<on_fast_path>(candidates, candidates_end, div_test::div_tests.data() + div_test::n_of_branchless_tests);
			count_passes(m += (candidates_end - candidates));



			// Do full primality tests, starting with base 2
			for (size_t* candidate_ptr = candidates; candidate_ptr < candidates_end; ++candidate_ptr)
			{
				const size_t candidate = *candidate_ptr;

				if (!franken::mpir_is_likely_prime_BPSW(candidate)) continue;
				count_passes(++n);

				// convert uint64_t to char array of ['0', '1'...] for MPIR
				char bin_str[64 + 1];
				auto result = std::to_chars(&bin_str[0], &bin_str[64], candidate, 2);
				*result.ptr = '\0';

				mpz_number.set_str(bin_str, 3);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;
				count_passes(++o);

				mpz_number.set_str(bin_str, 4);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;
				count_passes(++p);

				mpz_number.set_str(bin_str, 5);
				if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;
				count_passes(++q);

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

				mpz_number.set_str(bin_str, 13);
				if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 12); continue; }

				log_result(candidate, 13);
			}
		}

	private:
		gmp_randclass gmp_rand{ gmp_randinit_mt };
		mpz_class mpz_number = 0ull;

		count_passes(size_t a, ps15, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, passes, pc_hash = 0);
	};

} // namespace mbp



int main()
{
	mbp::print_preamble();

	mbp::find_multibase_primes mbp;
	mbp.run();
}
