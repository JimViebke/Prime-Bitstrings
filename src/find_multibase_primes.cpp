
#include <iomanip>
#include <iostream>
#include <sstream>

#include "find_multibase_primes.hpp"
#include "hardcoded_div_tests.hpp"
#include "io/io.hpp"
#include "merge_bitmasks.hpp"
#include "sieve.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/types.hpp"

namespace mbp
{
	static std::unique_ptr<std::array<sieve_container, prime_sieve::steps>> sieves =
		std::make_unique<decltype(sieves)::element_type>();
	static std::array<size_t, prime_sieve::steps> sieve_popcounts{};



	// 2x the expected number of candidates from the sieve passes
	constexpr size_t candidates_capacity = [] {
		double cleared = 0.0;
		for (size_t i = 1; i < small_primes_lookup.size(); ++i)
			cleared += (1.0 - cleared) * (1.0 / small_primes_lookup[i]);
		return 2 * size_t((1.0 - cleared) * sieve_container::size() * prime_sieve::steps);
	}();
	static std::array<uint64_t, candidates_capacity> candidates_storage alignas(64);



	// buffer candidates for full primality testing until we have 64
	constexpr size_t pt_buffer_capacity = 64;
	static std::array<uint64_t, pt_buffer_capacity> pt_buffer alignas(64);
	static size_t pt_buffer_size = 0;



	mbp::find_multibase_primes::find_multibase_primes()
	{
		gmp_rand.seed(mpir_ui{ 0xdeadbeef });

		count_passes(std::cout << "(counting passes)\n");
		count_passes(a = ps15 = b = c = d = e = bldt = bidt = 0);
		count_passes(b2 = b3 = b4 = b5 = passes = pc_hash = 0);
	}

	void mbp::find_multibase_primes::run()
	{
		constexpr size_t loop_size = 2ull * sieve_container::size() * prime_sieve::steps;

		size_t number = benchmark_mode ? bm_start : load_from_results();

		// Round starting number down to the nearest odd multiple of a product of primes
		number -= prime_sieve::product_of_static_sieve_primes; // n -= k
		number -= number % (2 * prime_sieve::product_of_static_sieve_primes); // n -= n % 2k
		number += prime_sieve::product_of_static_sieve_primes; // n += k

		// Align sieve on a multiple of 8, plus 1
		while ((number & 0b1111) != 0b0001)
			number -= (2 * prime_sieve::product_of_static_sieve_primes);

		prime_sieve::set_up_sieve_offsets_cache(number);

		size_t next_div_test_reorder = number + div_test::reorder_interval;

		// Start the clock after setup
		const auto start = util::current_time_in_ms();

		// (condition should optimize out)
		while (benchmark_mode ? number < bm_stop : true)
		{
			// Merge static sieve, popcount, gcd, and div test bitmasks
			for (size_t i = 0; i < prime_sieve::steps; ++i)
			{
				const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
				merge_bitmasks<1>(sieve_start, (*sieves)[i]);
			}
			for (size_t i = 0; i < prime_sieve::steps; ++i)
			{
				const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
				merge_bitmasks<2>(sieve_start, (*sieves)[i]);
			}
			for (size_t i = 0; i < prime_sieve::steps; ++i)
			{
				const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
				sieve_popcounts[i] = merge_bitmasks<3>(sieve_start, (*sieves)[i]);
				count_passes(a += sieve_popcounts[i]);
			}

			// Sieve until one of our density thresholds is reached
			for (size_t i = 0; i < prime_sieve::steps; ++i)
			{
				prime_sieve::partial_sieve((*sieves)[i], sieve_popcounts[i]);
				count_passes(ps15 += (*sieves)[i].count_bits());
			}

			uint64_t* const candidates = candidates_storage.data();
			uint64_t* candidates_end = candidates;

			// Convert 1-bit candidates to 64-bit candidates
			for (size_t i = 0; i < prime_sieve::steps; ++i)
			{
				const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
				candidates_end = prime_sieve::gather_sieve_results(candidates_end, (*sieves)[i], sieve_start);
			}

			// The upper 32 bits of a 64 bit integer only change every 4 billion ints.
			// Detect iterations where the upper bits can not change, and allow
			// functions to optimize based on this.
			if (util::upper_32_bits_match(number, number + loop_size))
			{
				candidates_end = div_tests<true>(candidates_end);
			}
			else
			{
				candidates_end = div_tests<false>(candidates_end);
			}

			// Collect prime candidates until we have filled our buffer, then do full primality tests starting with base 2
			for (const auto* ptr = candidates; ptr != candidates_end; ++ptr)
			{
				pt_buffer[pt_buffer_size++] = *ptr;

				if (pt_buffer_size == pt_buffer_capacity)
				{
					full_primality_tests(pt_buffer.data(), pt_buffer.data() + pt_buffer.size());
					pt_buffer_size = 0;
				}
			}



			number += loop_size;

			if (next_div_test_reorder <= number)
			{
				full_div_tests.update_div_test_order();
				next_div_test_reorder += div_test::reorder_interval;

			#if analyze_div_tests
				full_div_tests.print_div_tests();
			#endif
			}

			count_passes(++passes);
		} // end main loop



	#if analyze_div_tests
		full_div_tests.print_div_tests();
	#endif

		std::cout << "Finished. " << util::current_time_in_ms() - start << " ms elapsed\n";

		count_passes(std::cout << passes << " main loop iters\n");

		log_pass_counts("Passed bitmasks:       ", a, (bm_size / 2));
		log_pass_counts("Passed sieve:          ", ps15, a);
		log_pass_counts("Passed 4-rem tests:    ", b, ps15);
		log_pass_counts("Passed 8-rem tests:    ", c, b);
		log_pass_counts("Passed 12-rem tests:   ", d, c);
		log_pass_counts("Passed 16-rem tests:   ", e, d);
		log_pass_counts("P. branchless divtests:", bldt, e);
		log_pass_counts("P. branching divtests: ", bidt, bldt);
		log_pass_counts("Passed b2 BPSW test:   ", b2, bidt);
		log_pass_counts("Passed b3 prime test:  ", b3, b2);
		log_pass_counts("Passed b4 prime test:  ", b4, b3);
		log_pass_counts("Passed b5 prime test:  ", b5, b4);

		count_passes(std::cout << "\nhash of pass counts: " <<
					 std::hex << pc_hash << std::dec << '\n');
	}

	template<bool on_fast_path>
	uint64_t* mbp::find_multibase_primes::div_tests(uint64_t* candidates_end)
	{
		uint64_t* const candidates = candidates_storage.data();

		// Perform some div tests separately when a specialized implementation is faster

		// base 12 mod 29, and 6 mod 37 (4 remainders, part 2)
		candidates_end = two_div_tests_with_four_rems<12, 29, 6, 37, on_fast_path>(candidates, candidates_end);
		count_passes(b += (candidates_end - candidates));

		// bases 8 and 9 mod 17 (8 remainders)
		candidates_end = div_tests_with_8_rems<on_fast_path>(candidates, candidates_end);
		count_passes(c += (candidates_end - candidates));

		// bases 6, 7, and 11 mod 13 (12 remainders)
		candidates_end = div_tests_with_12_rems<on_fast_path>(candidates, candidates_end);
		count_passes(d += (candidates_end - candidates));

		// bases 3, 5, 6, 7, 10, 11 and 12 mod 17 (16 remainders)
		candidates_end = div_tests_with_16_rems<on_fast_path>(candidates, candidates_end);
		count_passes(e += (candidates_end - candidates));



		// Check for small prime factors across all bases
		candidates_end = full_div_tests.branchless_div_tests<on_fast_path>(candidates, candidates_end, div_test::n_of_branchless_tests);
		count_passes(bldt += (candidates_end - candidates));

		candidates_end = full_div_tests.branching_div_tests<on_fast_path>(candidates, candidates_end, div_test::n_of_branchless_tests);
		count_passes(bidt += (candidates_end - candidates));



		return candidates_end;
	}



	void print_config()
	{
		std::stringstream ss{};

		if constexpr (benchmark_mode)
		{
			ss << "Benchmarking from ";
			if constexpr (bm_start == p11) ss << "p11";
			else if constexpr (bm_start == p12) ss << "p12";
			else ss << bm_start;
			ss << ", size: " << bm_size / 1'000'000'000 << " B\n";
		}

		using namespace prime_sieve;
		ss << "Static sieve size: " << static_sieve_size
			<< ", primes: 3-" << static_sieve_primes.back() << '\n';
		ss << "Sieve limit: " << largest_sieve_prime
			<< ", vector/scalar thresholds: " << vector_density_threshold << ", " << scalar_density_threshold
			<< ", steps: " << steps
			<< ", candidate capacity: " << candidates_capacity << '\n';

		//ss << div_test::div_tests.size() << " div tests ("
		//	<< div_test::n_of_branchless_tests << "branchless + "
		//	<< div_test::div_tests.size() - div_test::n_of_branchless_tests << " branching), "
		//	<< div_test::n_of_primes << " prime factors, "
		//	<< "bases 2-" << div_test::up_to_base
		//	<< ", partial reorder every " << div_test::reorder_interval / 1'000'000'000 << " B\n";

		ss << "SPRP rounds: " << prime_test::n_random_bases << ", td limit: " << largest_sieve_prime << '\n';

	#define stringify(macro) #macro

		std::cout << ss.str() << std::endl;
	}
}

template uint64_t* mbp::find_multibase_primes::div_tests<true>(uint64_t*);
template uint64_t* mbp::find_multibase_primes::div_tests<false>(uint64_t*);
