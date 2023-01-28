
#include <sstream>

#include "bit_pattern_tests.hpp"
#include "find_multibase_primes.hpp"
#include "hardcoded_div_tests.hpp"
#include "io/io.hpp"
#include "sieve.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/simd.hpp"
#include "util/types.hpp"

namespace mbp
{
	static const sieve_container static_sieve = prime_sieve::generate_static_sieve();
	static sieve_container sieve;

	__forceinline size_t pc_lookup_idx(const size_t number)
	{
		constexpr size_t outer_48_bits_mask = ~(0xFFFFull << 1);

		return pop_count(number & outer_48_bits_mask) - 2; // -2 to normalize popcount 2-48 to idx 0-46
	}

	__forceinline size_t gcd_lookup_idx(const size_t number)
	{
		constexpr size_t outer_48_bits_mask = ~(0xFFFFull << 1);
		constexpr size_t even_mask = outer_48_bits_mask & 0x5555555555555555;
		constexpr size_t odd_mask = outer_48_bits_mask & 0xAAAAAAAAAAAAAAAA;

		const auto outer_even_pc = pop_count(number & even_mask);
		const auto outer_odd_pc = pop_count(number & odd_mask);
		return outer_even_pc - outer_odd_pc + 23; // +23 to normalize -23,24 to 0,47
	}

	template<size_t base, size_t prime>
	__forceinline size_t get_lookup_idx_for(const uint64_t number)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t outer_48_bits_mask = ~(0xFFFFull << 1);
		constexpr size_t max_sum_of_rems = get_sum_of_rems<prime, in_base<base>>(outer_48_bits_mask);

		const size_t sum = get_sum_of_rems<prime, in_base<base>>(number & outer_48_bits_mask);

		__assume(sum <= max_sum_of_rems);
		return sum % prime;
	}

	__forceinline void set_lookup_ptrs(const uint64_t next_number,
									   const uint8_t*& pc_lookup_ptr,
									   const uint8_t*& gcd_lookup_ptr,
									   const uint8_t*& b3m5_lookup_ptr,
									   const uint8_t*& b3m7_lookup_ptr,
									   const uint8_t*& b4m7_lookup_ptr,
									   const uint8_t*& b5m7_lookup_ptr)
	{
		constexpr uint64_t bits_1_16_mask = 0xFFFFull << 1;

		const size_t pc_idx = pc_lookup_idx(next_number);
		const size_t gcd_idx = gcd_lookup_idx(next_number);
		const size_t b3m5_idx = get_lookup_idx_for<3, 5>(next_number);
		const size_t b3m7_idx = get_lookup_idx_for<3, 7>(next_number);
		const size_t b4m7_idx = get_lookup_idx_for<4, 7>(next_number);
		const size_t b5m7_idx = get_lookup_idx_for<5, 7>(next_number);

		const uint64_t bits_1_16 = (next_number & bits_1_16_mask) >> 1;
		const size_t bit_offset = bits_1_16 & 0b111;
		const size_t byte_offset = bits_1_16 / 8;

		pc_lookup_ptr = pc_lookup[bit_offset][pc_idx].data() + byte_offset;
		gcd_lookup_ptr = gcd_lookup[bit_offset][gcd_idx].data() + byte_offset;
		b3m5_lookup_ptr = b3m5_lookup[bit_offset][b3m5_idx].data() + byte_offset;
		b3m7_lookup_ptr = b3m7_lookup[bit_offset][b3m7_idx].data() + byte_offset;
		b4m7_lookup_ptr = b4m7_lookup[bit_offset][b4m7_idx].data() + byte_offset;
		b5m7_lookup_ptr = b5m7_lookup[bit_offset][b5m7_idx].data() + byte_offset;
	}

	__forceinline void merge_one_block(uint64_t& number,
									   const uint8_t*& in,
									   uint8_t*& out,
									   const uint8_t*& pc_ptr,
									   const uint8_t*& gcd_ptr,
									   const uint8_t*& b3m5_ptr,
									   const uint8_t*& b3m7_ptr,
									   const uint8_t*& b4m7_ptr,
									   const uint8_t*& b5m7_ptr,
									   size_t& sieve_popcount)
	{
		constexpr uint64_t bits_1_16_mask = 0xFFFFull << 1;

		const size_t elements_to_rollover = (pow_2_16 - ((number & bits_1_16_mask) >> 1)) % pow_2_16; // map 65,536 -> 0

		uint64_t mask = *(uint64_t*)in;
		in += sizeof(uint64_t);

		const uint64_t bit_patterns_mask = (*(uint64_t*)pc_ptr &
											*(uint64_t*)gcd_ptr &
											*(uint64_t*)b3m5_ptr &
											*(uint64_t*)b3m7_ptr &
											*(uint64_t*)b4m7_ptr &
											*(uint64_t*)b5m7_ptr);
		pc_ptr += sizeof(uint64_t);
		gcd_ptr += sizeof(uint64_t);
		b3m5_ptr += sizeof(uint64_t);
		b3m7_ptr += sizeof(uint64_t);
		b4m7_ptr += sizeof(uint64_t);
		b5m7_ptr += sizeof(uint64_t);

		if (elements_to_rollover >= 64) // copy another block using lookup data
		{
			mask &= bit_patterns_mask;
		}
		else // elements_to_rollover is 0 to 63
		{
			set_lookup_ptrs(number + (elements_to_rollover * 2), pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr);
			const uint64_t new_mask = (*(uint64_t*)pc_ptr &
									   *(uint64_t*)gcd_ptr &
									   *(uint64_t*)b3m5_ptr &
									   *(uint64_t*)b3m7_ptr &
									   *(uint64_t*)b4m7_ptr &
									   *(uint64_t*)b5m7_ptr) << elements_to_rollover;

			const uint64_t select_from_old = (1ull << elements_to_rollover) - 1; // up to and including rollover
			const uint64_t select_from_new = ~select_from_old;

			mask &= ((bit_patterns_mask & select_from_old) |
					 (new_mask & select_from_new));

			// (re)set
			set_lookup_ptrs(number + (64ull * 2), pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr);
		}

		*(uint64_t*)out = mask;
		out += sizeof(uint64_t);

		sieve_popcount += pop_count(mask);

		number += (64ull * 2);
	}

	size_t copy_static_sieve_with_bit_pattern_filters(size_t number)
	{
		// 64 bits == 47 high bits + (16 "inner" bits) + 1 low bit (always set)
		constexpr size_t bits_1_16_mask = 0xFFFFull << 1;

		constexpr static uint256_t static_pc_shuf_lookup{ .m256i_u8{
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 } };
		constexpr static uint256_t static_nybble_mask{ .m256i_u64{
			0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F } };

		size_t sieve_popcount = 0;

		const uint8_t* in = static_sieve.data();
		const uint8_t* const aligned_end = in + ((sieve.size() / 256) * sizeof(uint256_t));

		uint8_t* out = sieve.data();

		// while in != aligned_end:
		// - iterate by 4 blocks until we hit aligned_end, or a rollover of bits 1-16, whichever happens first
		// - if we broke before aligned_end, we hit a rollover
		//   - handle the next 4 blocks seperately, then continue
		// handle cleanup

		// pointers into our large bitmask lookups
		const uint8_t* pc_ptr{};
		const uint8_t* gcd_ptr{};
		const uint8_t* b3m5_ptr{};
		const uint8_t* b3m7_ptr{};
		const uint8_t* b4m7_ptr{};
		const uint8_t* b5m7_ptr{};

		set_lookup_ptrs(number, pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr);

		for (; in != aligned_end; )
		{
			// The hot loop handles 256 elements (256 bits) per iteration.
			// How many times can we do this before reaching the lookups' ends?
			// Calculate this up front so the hot loop can run using a simpler condition.

			const size_t bits_1_16 = (number & bits_1_16_mask) >> 1;
			const size_t elements_to_rollover = (1ull << 16) - bits_1_16;

			constexpr size_t bytes_per_step = sizeof(uint256_t);
			constexpr size_t elements_per_step = bytes_per_step * 8;

			const size_t n_steps = std::min(elements_to_rollover / elements_per_step,
											(aligned_end - in) / bytes_per_step);

			const uint256_t nybble_mask = _mm256_loadu_si256(&static_nybble_mask);
			const uint256_t pc_shuf_lookup = _mm256_loadu_si256(&static_pc_shuf_lookup);
			uint256_t pc{};

			for (size_t offset = 0; offset < n_steps * bytes_per_step; offset += bytes_per_step)
			{
				const uint256_t pc_data = _mm256_loadu_si256((uint256_t*)(pc_ptr + offset));
				const uint256_t gcd_data = _mm256_loadu_si256((uint256_t*)(gcd_ptr + offset));
				const uint256_t b3m5_data = _mm256_loadu_si256((uint256_t*)(b3m5_ptr + offset));
				const uint256_t b3m7_data = _mm256_loadu_si256((uint256_t*)(b3m7_ptr + offset));
				const uint256_t b4m7_data = _mm256_loadu_si256((uint256_t*)(b4m7_ptr + offset));
				const uint256_t b5m7_data = _mm256_loadu_si256((uint256_t*)(b5m7_ptr + offset));
				const uint256_t ss_data = _mm256_loadu_si256((uint256_t*)(in + offset));

				uint256_t m1 = _mm256_and_si256(pc_data, gcd_data);
				uint256_t m2 = _mm256_and_si256(b3m5_data, b3m7_data);
				uint256_t m3 = _mm256_and_si256(b4m7_data, b5m7_data);
				uint256_t merged_data = _mm256_and_si256(_mm256_and_si256(m1, m2),
														 _mm256_and_si256(m3, ss_data));

				_mm256_storeu_si256((uint256_t*)(out + offset), merged_data);

				// data -> nybbles
				const uint256_t nybbles_lo = _mm256_and_si256(merged_data, nybble_mask);
				const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(merged_data, 4), nybble_mask);
				// nybbles -> 8-bit pcs
				uint256_t pc_256 = _mm256_add_epi8(_mm256_shuffle_epi8(pc_shuf_lookup, nybbles_lo),
												   _mm256_shuffle_epi8(pc_shuf_lookup, nybbles_hi));
				// 8-bit pcs -> 64-bit pcs
				pc_256 = _mm256_sad_epu8(pc_256, _mm256_setzero_si256());
				// 64-bit pcs -> running total
				pc = _mm256_add_epi64(pc, pc_256);
			}

			in += n_steps * bytes_per_step;
			out += n_steps * bytes_per_step;
			pc_ptr += n_steps * bytes_per_step;
			gcd_ptr += n_steps * bytes_per_step;
			b3m5_ptr += n_steps * bytes_per_step;
			b3m7_ptr += n_steps * bytes_per_step;
			b4m7_ptr += n_steps * bytes_per_step;
			b5m7_ptr += n_steps * bytes_per_step;

			number += n_steps * elements_per_step * 2; // * 2 because the sieve only contains odd numbers

			sieve_popcount += pc.m256i_u64[0];
			sieve_popcount += pc.m256i_u64[1];
			sieve_popcount += pc.m256i_u64[2];
			sieve_popcount += pc.m256i_u64[3];

			// if we are not at aligned_end, we are at a rollover
			if (in != aligned_end)
			{
				// handle 4 64-bit blocks
				for (size_t block = 0; block < 4; ++block)
				{
					merge_one_block(number, in, out, pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr, sieve_popcount);
				}
			}
		} // end main copy/merge loop

		// Less than 256 elements left. Generate instructions for 0-4 final copies.
		constexpr size_t leftover_elements = sieve.size() % 256;

		if constexpr (leftover_elements >= 64ull * 0)
		{
			merge_one_block(number, in, out, pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr, sieve_popcount);
		}
		if constexpr (leftover_elements >= 64ull * 1)
		{
			merge_one_block(number, in, out, pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr, sieve_popcount);
		}
		if constexpr (leftover_elements >= 64ull * 2)
		{
			merge_one_block(number, in, out, pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr, sieve_popcount);
		}
		if constexpr (leftover_elements >= 64ull * 3)
		{
			merge_one_block(number, in, out, pc_ptr, gcd_ptr, b3m5_ptr, b3m7_ptr, b4m7_ptr, b5m7_ptr, sieve_popcount);
		}

		return sieve_popcount;
	}



	// 2x the expected number of candidates from the sieve passes
	constexpr size_t candidates_capacity = [] {
		double cleared = 0.0;
		for (size_t i = 1; i < small_primes_lookup.size(); ++i)
			cleared += (1.0 - cleared) * (1.0 / small_primes_lookup[i]);
		return 2 * size_t((1.0 - cleared) * sieve.size() * prime_sieve::steps);
	}();
	static alignas(64) std::array<uint64_t, candidates_capacity> candidates_storage;

	mbp::find_multibase_primes::find_multibase_primes()
	{
		gmp_rand.seed(mpir_ui{ 0xdeadbeef });

		count_passes(std::cout << "(counting passes)\n");
		count_passes(a = ps15 = b = c = d = e = f = g = h = i = j = k = 0);
		count_passes(bldt = bidt = b2 = b3 = b4 = b5 = passes = pc_hash = 0);
	}

	void mbp::find_multibase_primes::run()
	{
		constexpr size_t loop_size = 2ull * sieve.size() * prime_sieve::steps;

		size_t number = benchmark_mode ? bm_start : load_from_results();

		// Round starting number down to the nearest odd multiple of the sieve size
		number -= sieve.size(); // n -= k
		number -= number % (2 * sieve.size()); // n -= n % 2k
		number += sieve.size(); // n += k

		prime_sieve::set_up_sieve_offsets_cache(number);

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
				div_test::update_div_test_order();
				next_div_test_reorder += div_test::reorder_interval;

			#if analyze_div_tests
				div_test::print_div_tests();
				//run_div_test_analysis(number);
			#endif
			}

			count_passes(++passes);
		} // end main loop



	#if analyze_div_tests
		  // Run one final step before exiting
		  // run_div_test_analysis();
		div_test::print_div_tests();
	#endif

		std::cout << "Finished. " << util::current_time_in_ms() - start << " ms elapsed\n";

		count_passes(std::cout << passes << " main loop iters\n");

		log_pass_counts("Passed static sieve and\n"\
						"  bit pattern filters: ", a, (bm_size / 2));
		log_pass_counts("Passed sieve:          ", ps15, a);
		log_pass_counts("Passed 4-rem tests:    ", b, ps15);
		log_pass_counts("Passed b13m17 test:    ", c, b);
		log_pass_counts("Passed 4-rem tests, p2:", d, c);
		log_pass_counts("Passed 5-rem tests:    ", e, d);
		log_pass_counts("Passed 6-rem tests:    ", f, e);
		log_pass_counts("Passed 10-rem tests:   ", g, f);
		log_pass_counts("Passed 3-rem tests:    ", h, g);
		log_pass_counts("Passed 8-rem tests:    ", i, h);
		log_pass_counts("Passed 12-rem tests:   ", j, i);
		log_pass_counts("Passed 16-rem tests:   ", k, j);
		log_pass_counts("P. branchless divtests:", bldt, k);
		log_pass_counts("P. branching divtests: ", bidt, bldt);
		log_pass_counts("Passed b2 BPSW test:   ", b2, bidt);
		log_pass_counts("Passed b3 prime test:  ", b3, b2);
		log_pass_counts("Passed b4 prime test:  ", b4, b3);
		log_pass_counts("Passed b5 prime test:  ", b5, b4);

		count_passes(std::cout << "\nhash of pass counts: " <<
					 std::hex << pc_hash << std::dec << '\n');
	}

	template<bool on_fast_path>
	void mbp::find_multibase_primes::main_loop(const uint64_t number)
	{
		uint64_t* const candidates = candidates_storage.data();
		uint64_t* candidates_end = candidates;

		for (size_t sieve_step = 0; sieve_step < prime_sieve::steps; ++sieve_step)
		{
			const uint64_t sieve_start = number + (sieve_step * sieve.size() * 2);

			// Merge the static sieve, popcount, and gcd data
			const size_t sieve_popcount = copy_static_sieve_with_bit_pattern_filters(sieve_start);
			count_passes(a += sieve_popcount);

			prime_sieve::partial_sieve(sieve_start, sieve, sieve_popcount);
			count_passes(ps15 += sieve.count_bits());

			candidates_end = prime_sieve::gather_sieve_results(
				candidates_end, sieve, sieve_start);
		}



		// Perform some div tests separately when a specialized implementation is faster

		// base 4 mod 17, and bases 5 and 8 mod 13 (4 remainders)
		candidates_end = div_tests_with_four_rems<on_fast_path>(candidates, candidates_end);
		count_passes(b += (candidates_end - candidates));

		// base 13 mod 17 (4 remainders, 8 candidates)
		candidates_end = base13_mod17_div_test<on_fast_path>(candidates, candidates_end);
		count_passes(c += (candidates_end - candidates));

		// slower:
		// candidates_end = two_div_tests_with_four_rems<5, 13, 8, 13, on_fast_path>(candidates, candidates_end);
		// count_passes(b += (candidates_end - candidates));
		// candidates_end = two_div_tests_with_four_rems<4, 17, 13, 17, on_fast_path>(candidates, candidates_end);
		// count_passes(c += (candidates_end - candidates));

		// also slower:
		// candidates_end = four_div_tests_with_four_rems<on_fast_path>(candidates, candidates_end);
		// count_passes(b += (candidates_end - candidates));

		// base 12 mod 29, and 6 mod 37 (4 remainders, part 2)
		candidates_end = two_div_tests_with_four_rems<12, 29, 6, 37, on_fast_path>(candidates, candidates_end);
		count_passes(d += (candidates_end - candidates));

		// bases 3, 4, 5, and 9 mod 11 (5 remainders)
		candidates_end = div_tests_with_five_rems<on_fast_path>(candidates, candidates_end);
		count_passes(e += (candidates_end - candidates));

		// bases 4 and 10 mod 13 (6 remainders)
		candidates_end = div_tests_with_six_rems<on_fast_path>(candidates, candidates_end);
		count_passes(f += (candidates_end - candidates));

		// bases 6, 7, and 8 mod 11 (10 remainders)
		candidates_end = div_tests_with_10_rems<on_fast_path>(candidates, candidates_end);
		count_passes(g += (candidates_end - candidates));

		// bases 3 and 9 mod 13 (3 remainders)
		candidates_end = div_tests_with_three_rems<on_fast_path>(candidates, candidates_end);
		count_passes(h += (candidates_end - candidates));

		// bases 8 and 9 mod 17 (8 remainders)
		candidates_end = div_tests_with_8_rems<on_fast_path>(candidates, candidates_end);
		count_passes(i += (candidates_end - candidates));

		// bases 6, 7, and 11 mod 13 (12 remainders)
		candidates_end = div_tests_with_12_rems<on_fast_path>(candidates, candidates_end);
		count_passes(j += (candidates_end - candidates));

		// bases 3, 5, 6, 7, 10, 11 and 12 mod 17 (16 remainders)
		candidates_end = div_tests_with_16_rems<on_fast_path>(candidates, candidates_end);
		count_passes(k += (candidates_end - candidates));

		// bases 8 and 12 mod 19 (6 remainders)
		//candidates_end = two_div_tests_with_six_rems<on_fast_path>(candidates, candidates_end);
		//count_passes(dt6b += (candidates_end - candidates));

		// bases 4, 5, 6, and 9 mod 19 (9 remainders)
		//candidates_end = div_tests_with_nine_rems<on_fast_path>(candidates, candidates_end);
		//count_passes(dt9 += (candidates_end - candidates));

		// bases 7 and 11 mod 19 (3 remainders)
		//candidates_end = two_div_tests_with_three_rems<on_fast_path>(candidates, candidates_end);
		//count_passes(dt3b += (candidates_end - candidates));



		// Check for small prime factors across all bases
		candidates_end = div_test::branchless_div_tests<on_fast_path>(candidates, candidates_end, div_test::n_of_branchless_tests);
		count_passes(bldt += (candidates_end - candidates));

		candidates_end = div_test::branching_div_tests<on_fast_path>(candidates, candidates_end, div_test::n_of_branchless_tests);
		count_passes(bidt += (candidates_end - candidates));



		// Do full primality tests, starting with base 2
		full_primality_tests(candidates, candidates_end);
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
		//	<< "bases 3-" << div_test::up_to_base
		//	<< ", partial reorder every " << div_test::reorder_interval / 1'000'000'000 << " B\n";

		ss << "SPRP rounds: " << prime_test::n_random_bases << ", td limit: " << largest_sieve_prime << '\n';

	#define stringify(macro) #macro

		std::cout << ss.str() << std::endl;
	}
}

template void mbp::find_multibase_primes::main_loop<true>(const size_t);
template void mbp::find_multibase_primes::main_loop<false>(const size_t);
