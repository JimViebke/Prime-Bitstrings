
#include <algorithm>
#include <iostream>

#include "../util/simd.hpp"
#include "multibase_div_tests.hpp"

#if analyze_div_tests
#include <iomanip>
#include <sstream>
#endif

namespace mbp::div_test
{
	namespace detail
	{
		std::vector<div_test_t> generate_div_tests_impl()
		{
			std::vector<uncompressed_div_test_t> uncompressed_dts;

			for (size_t i = 1; i < n_of_primes; ++i) // for each small prime starting from 3
			{
				for (size_t base = 2; base <= up_to_base; ++base) // for each base 2..n
				{
					const auto p = small_primes_lookup[i];

					// Skip always-composite cases
					if (base % p == 0) continue;

					// Skip candidates removed by the static sieve
					if (base == 2 && p <= prime_sieve::static_sieve_primes.back()) continue;

					// Implemented with bit masks
					// pass 1
					if (base == 3 && p == 5) continue;
					if (base == 4 && p == 17) continue;
					if (base == 5 && p == 13) continue;
					if (base == 8 && p == 13) continue;
					if (base == 13 && p == 17) continue;
					// pass 2
					if (base == 3 && p == 7) continue;
					if (base == 3 && p == 13) continue;
					if (base == 4 && p == 7) continue;
					if (base == 4 && p == 13) continue;
					if (base == 5 && p == 7) continue;
					if (base == 9 && p == 13) continue;
					if (base == 10 && p == 13) continue;
					// pass 3
					if (base == 3 && p == 11) continue;
					if (base == 4 && p == 11) continue;
					if (base == 5 && p == 11) continue;
					if (base == 9 && p == 11) continue;

					// Always suppress hardcoded div tests
					// 4 remainders
					if (base == 12 && p == 29) continue;
					if (base == 6 && p == 37) continue;

					// 10 remainders
					if (base == 6 && p == 11) continue;
					if (base == 7 && p == 11) continue;
					if (base == 8 && p == 11) continue;

					// 8 remainders
					if (base == 8 && p == 17) continue;
					if (base == 9 && p == 17) continue;

					// 12 remainders
					if (base == 6 && p == 13) continue;
					if (base == 7 && p == 13) continue;
					if (base == 11 && p == 13) continue;

					// 16 remainders
					if (base == 3 && p == 17) continue;
					if (base == 5 && p == 17) continue;
					if (base == 6 && p == 17) continue;
					if (base == 7 && p == 17) continue;
					if (base == 10 && p == 17) continue;
					if (base == 11 && p == 17) continue;
					if (base == 12 && p == 17) continue;

					// second set with 6 remainders
					//if (base == 8 && p == 19) continue;
					//if (base == 12 && p == 19) continue;

					// 9 remainders
					//if (base == 4 && p == 19) continue;
					//if (base == 5 && p == 19) continue;
					//if (base == 6 && p == 19) continue;
					//if (base == 9 && p == 19) continue;

					// second set with 3 remainders
					//if (base == 7 && p == 19) continue;
					//if (base == 11 && p == 19) continue;

				#if !analyze_div_tests or suppress_extra_div_tests
					if (base == 4 && p == 3) continue; //  base  4^n % 3 unused
					if (base == 5 && p == 3) continue; //  base  5^n % 3 unused
					if (base == 7 && p == 3) continue; //  base  7^n % 3 unused
					if (base == 8 && p == 3) continue; //  base  8^n % 3 unused
					if (base == 10 && p == 3) continue; // base 10^n % 3 unused
					if (base == 11 && p == 3) continue; // base 11^n % 3 unused
					if (base == 13 && p == 3) continue; // base 13^n % 3 unused

					if (base == 4 && p == 5) continue; //  base  4^n % 5 unused
					if (base == 6 && p == 5) continue; //  base  6^n % 5 unused
					if (base == 7 && p == 5) continue; //  base  7^n % 5 unused
					if (base == 9 && p == 5) continue; //  base  9^n % 5 unused
					if (base == 11 && p == 5) continue; // base 11^n % 5 unused
					if (base == 12 && p == 5) continue; // base 12^n % 5 unused
					if (base == 13 && p == 5) continue; // base 13^n % 5 unused

					if (base == 6 && p == 7) continue; //  base  6^n % 7 unused
					if (base == 8 && p == 7) continue; //  base  8^n % 7 unused
					if (base == 9 && p == 7) continue; //  base  9^n % 7 unused
					if (base == 13 && p == 7) continue; // base 13^n % 7 unused

					if (base == 10 && p == 11) continue; // base 10^n % 11 unused
					if (base == 12 && p == 11) continue; // base 12^n % 11 unused
					if (base == 13 && p == 11) continue; // base 13^n % 11 unused

					if (base == 12 && p == 13) continue; // base 12^n % 13 unused

					// If two div tests are effectively identical, remove one
					if (base == 8 && p == 5) continue; //  base  8^n % 5 is congruent to 3^n % 5
					if (base == 10 && p == 7) continue; // base 10^n % 7 is congruent to 3^n % 7
					if (base == 11 && p == 7) continue; // base 11^n % 7 is congruent to 4^n % 7
					if (base == 12 && p == 7) continue; // base 12^n % 7 is congruent to 5^n % 7

				#endif

					uncompressed_div_test_t dt{ .base = base_t(base), .prime_idx = prime_idx_t(i) };

					// calculate base^j mod prime, where j is the place value
					for (size_t j = 0; j < 64; ++j)
					{
						remainder_t rem = remainder_t(pk::powMod(base, j, p));
						if (rem == 1 && j > 0)
						{
							// The pattern is repeating; stop generating further terms
							break;
						}

						dt.remainders[j] = rem;
						dt.n_of_remainders++;
					}

					uncompressed_dts.push_back(dt);
				}
			}

			// move base 2 tests to the end
			std::stable_partition(uncompressed_dts.begin(),
								  uncompressed_dts.end(),
								  [](auto dt) { return dt.base != 2; });

			std::vector<div_test_t> div_tests;
			div_tests.reserve(uncompressed_dts.size());

			for (const auto& udt : uncompressed_dts)
			{
				div_test_t dt{
					.prime_idx = udt.prime_idx,
					.n_of_remainders = udt.n_of_remainders,
					.remainders = udt.remainders };

				// Repeat the remainders so all div tests have 64 terms
				for (size_t i = dt.n_of_remainders; i < 64; ++i)
				{
					dt.remainders[i] = dt.remainders[i - dt.n_of_remainders];
				}

			#if analyze_div_tests
				// We do need to copy this if we're in analyze mode
				dt.base = udt.base;
			#endif

				div_tests.push_back(dt);
			}

			return div_tests;
		}

		size_t calculate_prime_factor_lookup_size()
		{
			const div_tests_t div_tests_temp = detail::generate_div_tests_impl();

			size_t largest_sum = 0;

			// Calculate the largest possible sum of remainders of a number with every bit set
			for (const auto& div_test : div_tests_temp)
			{
				const size_t sum = std::accumulate(div_test.remainders.begin(), div_test.remainders.end(), size_t(0));

				if (sum > largest_sum)
					largest_sum = sum;
			}

			// + 1 so "lookup[sum]" is always in range
			return largest_sum + 1;
		}

		// Faster version
		std::vector<prime_lookup_t> build_prime_factor_lookup_old()
		{
			const size_t prime_factor_lookup_size = calculate_prime_factor_lookup_size();

			std::vector<prime_lookup_t> lookup;
			lookup.reserve(prime_factor_lookup_size);

			for (size_t i = 0; i < prime_factor_lookup_size; ++i)
			{
				prime_lookup_t entry = 0;
				for (prime_lookup_t p = 0; p < div_test::n_of_primes; ++p)
				{
					entry |= (prime_lookup_t((i % small_primes_lookup[p]) == 0) << p);
				}

				lookup.push_back(entry);
			}

			return lookup;
		}

		// Slower version
		std::vector<prime_lookup_t> build_prime_factor_lookup_new()
		{
			const size_t prime_factor_lookup_size = calculate_prime_factor_lookup_size();
			std::vector<prime_lookup_t> lookup(prime_factor_lookup_size, 0);

			for (size_t i = 0; i < div_test::n_of_primes; ++i)
			{
				const prime_lookup_t p = small_primes_lookup[i];

				for (prime_lookup_t j = 0; j < prime_factor_lookup_size; j += p)
				{
					lookup[j] |= (1ull << i);
				}
			}

			return lookup;
		}

		std::array<std::vector<uint8_t>, n_of_primes> build_indivisible_lookup()
		{
			const div_tests_t div_tests_temp = detail::generate_div_tests_impl();

			std::array<std::vector<uint8_t>, n_of_primes> lookup;

			// start from 5, at idx 2
			for (size_t prime_idx = 2; prime_idx < n_of_primes; ++prime_idx)
			{
				// calculate the largest possible sum of remainders, for any
				// div test using this prime

				size_t largest_sum = 0;

				for (const auto& div_test : div_tests_temp)
				{
					if (div_test.prime_idx != prime_idx) continue;

					const size_t sum = std::accumulate(div_test.remainders.begin(), div_test.remainders.end(), size_t(0));

					if (sum > largest_sum)
						largest_sum = sum;
				}

				if (largest_sum == 0)
				{
					// no div test uses this prime, generate remainders, then the sum of remainders, for each base

					const auto prime = small_primes_lookup[prime_idx];

					for (size_t base = 2; base < up_to_base; ++base)
					{
						std::array<remainder_t, 64> rems{};
						size_t i = 0;
						for (; i < 64; ++i)
						{
							remainder_t rem = remainder_t(pk::powMod(base, i, prime));
							// break when/if the pattern repeats
							if (rem == 1 && i > 0) break;
							rems[i] = rem; // save
						}

						// extend rems to 64
						for (size_t j = i; j < 64; ++j)
						{
							rems[j] = rems[j - i];
						}

						const size_t sum = std::accumulate(rems.begin(), rems.end(), size_t{ 0 });
						largest_sum = std::max(sum, largest_sum);
					}
				}

				std::vector<uint8_t>& v = lookup[prime_idx];
				v.assign(largest_sum + 1, sizeof(uint64_t));

				for (size_t i = 0; i < v.size(); i += small_primes_lookup[prime_idx])
				{
					v[i] = 0;
				}
			}

			return lookup;
		}

		// increment_in_bytes = indivisible_by[prime_index][sum_of_remainders]
		const std::array<std::vector<uint8_t>, n_of_primes> indivisible_by = build_indivisible_lookup();
	}



	full_div_tests::full_div_tests() :
		div_tests{ detail::generate_div_tests_impl() },
		permuted_div_tests{}
	{
		const size_t n_of_branching_tests = div_tests.size() - n_of_branchless_tests;

		if (n_of_branching_tests & 1)
		{
			std::cout << "error - odd number of branching div tests";
			std::cin.ignore();
		}

		permuted_div_tests.assign(n_of_branching_tests, decltype(permuted_div_tests)::value_type{});
		permute_div_tests();
	}



	template<bool on_fast_path>
	uint64_t* full_div_tests::branchless_div_tests(uint64_t* const candidates_begin,
												   uint64_t* const candidates_end,
												   const size_t n_of_tests)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr uint64_t upper_bits_mask = uint64_t(-1) << 32;

		uint64_t* shrinking_end = candidates_end;

		static constexpr uint64_t static_shuffle_mask_byte_0[4] = {
			0x0000000000000000, 0x0808080808080808, 0x0000000000000000, 0x0808080808080808 };
		static constexpr uint64_t static_shuffle_mask_byte_1[4] = {
			0x0101010101010101, 0x0909090909090909, 0x0101010101010101, 0x0909090909090909 };
		static constexpr uint64_t static_shuffle_mask_byte_2[4] = {
			0x0202020202020202, 0x0A0A0A0A0A0A0A0A, 0x0202020202020202, 0x0A0A0A0A0A0A0A0A };
		static constexpr uint64_t static_shuffle_mask_byte_3[4] = {
			0x0303030303030303, 0x0B0B0B0B0B0B0B0B, 0x0303030303030303, 0x0B0B0B0B0B0B0B0B };

		const uint256_t shuffle_mask_byte_0 = _mm256_load_si256((uint256_t*)static_shuffle_mask_byte_0);
		const uint256_t shuffle_mask_byte_1 = _mm256_load_si256((uint256_t*)static_shuffle_mask_byte_1);
		const uint256_t shuffle_mask_byte_2 = _mm256_load_si256((uint256_t*)static_shuffle_mask_byte_2);
		const uint256_t shuffle_mask_byte_3 = _mm256_load_si256((uint256_t*)static_shuffle_mask_byte_3);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		for (size_t i = 0; i < n_of_tests; ++i)
		{
			div_test_t& div_test = div_tests[i];

			const uint64_t* input = candidates_begin;
			uint64_t* output = candidates_begin;

			uint64_t upper_bits = *input & upper_bits_mask;
			const uint8_t* indivisible_ptr = indivisible_by[div_test.prime_idx].data() +
				util::vcl_hadd_x(_mm256_and_si256(util::expand_bits_to_bytes(upper_bits >> 32),
												  _mm256_loadu_si256((uint256_t*)&div_test.remainders[32])));

			if constexpr (on_fast_path)
			{
				const uint256_t rems_0 = _mm256_set1_epi64x(*(uint64_t*)&div_test.remainders[8 * 0]);
				const uint256_t rems_1 = _mm256_set1_epi64x(*(uint64_t*)&div_test.remainders[8 * 1]);
				const uint256_t rems_2 = _mm256_set1_epi64x(*(uint64_t*)&div_test.remainders[8 * 2]);
				const uint256_t rems_3 = _mm256_set1_epi64x(*(uint64_t*)&div_test.remainders[8 * 3]);

				// calculate the number of candidates, rounded down to the nearest 4
				const size_t n_of_candidates = shrinking_end - input;
				const uint64_t* const rounded_end = input + (n_of_candidates - (n_of_candidates % 4));

				// load four candidates
				uint256_t candidates = _mm256_loadu_si256((uint256_t*)input);

				alignas(32) volatile uint64_t sums[4]{};

				// run vector instructions one iteration ahead
				{
					// convert bits to bytes, storing 8 copies of the nth byte of four candidates per register
					uint256_t candidates_byte_0 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_0);
					uint256_t candidates_byte_1 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_1);
					uint256_t candidates_byte_2 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_2);
					uint256_t candidates_byte_3 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_3);

					// load ahead
					candidates = _mm256_loadu_si256((uint256_t*)(input + 4));

					// convert bits to bytes, cont'd
					candidates_byte_0 = _mm256_andnot_si256(candidates_byte_0, and_mask);
					candidates_byte_1 = _mm256_andnot_si256(candidates_byte_1, and_mask);
					candidates_byte_2 = _mm256_andnot_si256(candidates_byte_2, and_mask);
					candidates_byte_3 = _mm256_andnot_si256(candidates_byte_3, and_mask);
					candidates_byte_0 = _mm256_cmpeq_epi8(candidates_byte_0, _mm256_setzero_si256());
					candidates_byte_1 = _mm256_cmpeq_epi8(candidates_byte_1, _mm256_setzero_si256());
					candidates_byte_2 = _mm256_cmpeq_epi8(candidates_byte_2, _mm256_setzero_si256());
					candidates_byte_3 = _mm256_cmpeq_epi8(candidates_byte_3, _mm256_setzero_si256());

					// mask to select remainders
					candidates_byte_0 = _mm256_and_si256(candidates_byte_0, rems_0);
					candidates_byte_1 = _mm256_and_si256(candidates_byte_1, rems_1);
					candidates_byte_2 = _mm256_and_si256(candidates_byte_2, rems_2);
					candidates_byte_3 = _mm256_and_si256(candidates_byte_3, rems_3);

					// vertically sum selected remainders
					uint256_t sums_01 = _mm256_add_epi8(candidates_byte_0, candidates_byte_1);
					uint256_t sums_23 = _mm256_add_epi8(candidates_byte_2, candidates_byte_3);
					// h-sum 4 sets of 8 consecutive bytes
					sums_01 = _mm256_sad_epu8(sums_01, _mm256_setzero_si256());
					sums_23 = _mm256_sad_epu8(sums_23, _mm256_setzero_si256());
					// final v-sum
					const uint256_t sums_0123 = _mm256_add_epi16(sums_01, sums_23);

					// store results to the stack
					_mm256_store_si256((uint256_t*)sums, sums_0123);
				}

				for (; input < rounded_end; )
				{
					// run vector instructions one iteration ahead
					uint256_t candidates_byte_0 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_0);
					uint256_t candidates_byte_1 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_1);
					uint256_t candidates_byte_2 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_2);
					uint256_t candidates_byte_3 = _mm256_shuffle_epi8(candidates, shuffle_mask_byte_3);

					// load two iterations ahead
					candidates = _mm256_loadu_si256((uint256_t*)(input + 8));

					candidates_byte_0 = _mm256_andnot_si256(candidates_byte_0, and_mask);
					candidates_byte_1 = _mm256_andnot_si256(candidates_byte_1, and_mask);
					candidates_byte_2 = _mm256_andnot_si256(candidates_byte_2, and_mask);
					candidates_byte_3 = _mm256_andnot_si256(candidates_byte_3, and_mask);
					candidates_byte_0 = _mm256_cmpeq_epi8(candidates_byte_0, _mm256_setzero_si256());
					candidates_byte_1 = _mm256_cmpeq_epi8(candidates_byte_1, _mm256_setzero_si256());
					candidates_byte_2 = _mm256_cmpeq_epi8(candidates_byte_2, _mm256_setzero_si256());
					candidates_byte_3 = _mm256_cmpeq_epi8(candidates_byte_3, _mm256_setzero_si256());

					candidates_byte_0 = _mm256_and_si256(candidates_byte_0, rems_0);
					candidates_byte_1 = _mm256_and_si256(candidates_byte_1, rems_1);
					candidates_byte_2 = _mm256_and_si256(candidates_byte_2, rems_2);
					candidates_byte_3 = _mm256_and_si256(candidates_byte_3, rems_3);

					uint256_t sums_01 = _mm256_add_epi8(candidates_byte_0, candidates_byte_1);
					uint256_t sums_23 = _mm256_add_epi8(candidates_byte_2, candidates_byte_3);
					sums_01 = _mm256_sad_epu8(sums_01, _mm256_setzero_si256());
					sums_23 = _mm256_sad_epu8(sums_23, _mm256_setzero_si256());
					const uint256_t sums_0123 = _mm256_add_epi16(sums_01, sums_23);

					const uint64_t inc_a = indivisible_ptr[sums[0]];
					const uint64_t inc_b = indivisible_ptr[sums[1]];
					const uint64_t inc_c = indivisible_ptr[sums[2]];
					const uint64_t inc_d = indivisible_ptr[sums[3]];

					// load and increment
					const uint64_t candidate_a = *input;
					const uint64_t candidate_b = *(input + 1);
					const uint64_t candidate_c = *(input + 2);
					const uint64_t candidate_d = *(input + 3);
					input += 4;

					*output = candidate_a; // always write
					// branchless conditional increment
					output = (uint64_t*)(((uint8_t*)output) + inc_a);

					*output = candidate_b;
					output = (uint64_t*)(((uint8_t*)output) + inc_b);

					*output = candidate_c;
					output = (uint64_t*)(((uint8_t*)output) + inc_c);

					*output = candidate_d;
					output = (uint64_t*)(((uint8_t*)output) + inc_d);

					// store the above results to the stack
					_mm256_store_si256((uint256_t*)sums, sums_0123);
				}

				// if there are remaining candidates, the loop below runs 1-3 times and handles them

			} // end if (on_fast_path)

			const uint256_t rems_lo = _mm256_loadu_si256((uint256_t*)&div_test.remainders[0]);
			const uint256_t rems_hi = _mm256_loadu_si256((uint256_t*)&div_test.remainders[32]);

			uint64_t number = *input; // load one iteration ahead

			for (; input < shrinking_end; )
			{
				if constexpr (!on_fast_path)
				{
					// Upper 32 bits only change every 4B integers - re-calculate when stale
					if ((number & upper_bits_mask) != upper_bits)
					{
						upper_bits = number & upper_bits_mask;
						indivisible_ptr = indivisible_by[div_test.prime_idx].data() +
							util::vcl_hadd_x(_mm256_and_si256(util::expand_bits_to_bytes(upper_bits >> 32), rems_hi));
					}
				}

				const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));

				*output = number; // always write
				number = *++input; // load one iteration ahead

				const uint256_t rems = _mm256_and_si256(mask_lower, rems_lo);
				const uint64_t rem = util::vcl_hadd_x(rems);

				// increment only if we haven't found a prime factor yet
				const auto inc = indivisible_ptr[rem];
				output = (uint64_t*)(((uint8_t*)output) + inc);
			}

			div_test.hits += uint32_t(input - output);

			// output is one past the last element; use it as the new end
			shrinking_end = output;
		}

		return shrinking_end;
	}

	template uint64_t* full_div_tests::branchless_div_tests<true>(uint64_t* const, uint64_t* const, const size_t);
	template uint64_t* full_div_tests::branchless_div_tests<false>(uint64_t* const, uint64_t* const, const size_t);



	template<bool on_fast_path>
	uint64_t* full_div_tests::branching_div_tests(uint64_t* input,
												  const uint64_t* const candidates_end,
												  const size_t start_offset)
	{
		using namespace div_test;

		static_assert(sizeof(remainder_t) == 1);

		constexpr uint64_t upper_bits_mask = uint64_t(-1) << 32;

		div_test_t* div_tests_start = div_tests.data() + start_offset;
		const div_test_t* const div_tests_end = div_tests.data() + div_tests.size();

		// store the bitstring across four registers, where high and low lanes store the same 16 bytes
		uint64_t upper_bits = *input & upper_bits_mask;
		uint256_t mask_upper = util::expand_bits_to_bytes(upper_bits >> 32);
		uint256_t candidate_bytes_45 = _mm256_permute4x64_epi64(mask_upper, 0b01'00'01'00); // 1, 0, 1, 0
		uint256_t candidate_bytes_67 = _mm256_permute4x64_epi64(mask_upper, 0b11'10'11'10); // 3, 2, 3, 2

		uint64_t* output = input;

		for (; input < candidates_end; ++input)
		{
			const uint64_t number = *input;
			*output = number; // always write

			if constexpr (!on_fast_path)
			{
				// Upper 32 bits only change every 4B integers - re-calculate when stale
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;
					mask_upper = util::expand_bits_to_bytes(upper_bits >> 32);
					candidate_bytes_45 = _mm256_permute4x64_epi64(mask_upper, 0b01'00'01'00); // 1, 0, 1, 0
					candidate_bytes_67 = _mm256_permute4x64_epi64(mask_upper, 0b11'10'11'10); // 3, 2, 3, 2
				}
			}

			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));
			const uint256_t candidate_bytes_01 = _mm256_permute4x64_epi64(mask_lower, 0b01'00'01'00); // 1, 0, 1, 0
			const uint256_t candidate_bytes_23 = _mm256_permute4x64_epi64(mask_lower, 0b11'10'11'10); // 3, 2, 3, 2

			size_t is_candidate = 1;

			// load 2x 64 remainders one iteration ahead
			const auto* permuted_ptr = permuted_div_tests.data();
			uint256_t rems_0 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 0)));
			uint256_t rems_1 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 1)));
			uint256_t rems_2 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 2)));
			uint256_t rems_3 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 3)));

			for (div_test_t* div_test_ptr = div_tests_start; div_test_ptr < div_tests_end; )
			{
				div_test_t& div_test_0 = *(div_test_ptr + 0);
				div_test_t& div_test_1 = *(div_test_ptr + 1);

				// mask 64 digits against 2 div tests
				uint256_t selected_rems_0 = _mm256_and_si256(rems_0, candidate_bytes_01);
				uint256_t selected_rems_1 = _mm256_and_si256(rems_1, candidate_bytes_23);
				uint256_t selected_rems_2 = _mm256_and_si256(rems_2, candidate_bytes_45);
				uint256_t selected_rems_3 = _mm256_and_si256(rems_3, candidate_bytes_67);

				// load 2x 64 remainders one iteration ahead
				permuted_ptr += 2;
				div_test_ptr += 2;
				rems_0 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 0)));
				rems_1 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 1)));
				rems_2 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 2)));
				rems_3 = _mm256_loadu_si256((uint256_t*)(permuted_ptr->data() + (32 * 3)));

				// horizontally add sets of 8 consecutive remainders
				selected_rems_0 = _mm256_sad_epu8(selected_rems_0, _mm256_setzero_si256());
				selected_rems_1 = _mm256_sad_epu8(selected_rems_1, _mm256_setzero_si256());
				selected_rems_2 = _mm256_sad_epu8(selected_rems_2, _mm256_setzero_si256());
				selected_rems_3 = _mm256_sad_epu8(selected_rems_3, _mm256_setzero_si256());

				// vertically add sets of 16-bit remainders
				uint256_t sum = _mm256_add_epi64(
					_mm256_add_epi64(selected_rems_0, selected_rems_1),
					_mm256_add_epi64(selected_rems_2, selected_rems_3));

				// final horizontal add
				sum = _mm256_add_epi64(sum, _mm256_srli_si256(sum, 8));

				const uint64_t sum_0 = _mm256_extract_epi64(sum, 0);

				if (has_small_prime_factor(sum_0, div_test_0.prime_idx))
				{
					div_test_0.hits++;
					is_candidate = 0;
					break;
				}

				const uint64_t sum_1 = _mm256_extract_epi64(sum, 2); // extract from bottom of high lane
				if (has_small_prime_factor(sum_1, div_test_1.prime_idx))
				{
					div_test_1.hits++;
					is_candidate = 0;
					break;
				}
			}

			// conditionally increment
			output += is_candidate;
		}

		return output;
	}

	template uint64_t* full_div_tests::branching_div_tests<true>(uint64_t*, const uint64_t* const, const size_t);
	template uint64_t* full_div_tests::branching_div_tests<false>(uint64_t*, const uint64_t* const, const size_t);



	void full_div_tests::permute_div_tests()
	{
		// (re)build a list of div tests, where pairs of div tests are permuted and stored across four registers
		auto* permuted_dt = permuted_div_tests.data();
		for (auto it = div_tests.cbegin() + n_of_branchless_tests; it < div_tests.cend(); it += 2, permuted_dt += 2)
		{
			// _mm256_loadu2_m128i() takes its args in (high, low) order
			const uint256_t rems_0 = _mm256_loadu2_m128i((uint128_t*)&(it + 1)->remainders[16 * 0],
														 (uint128_t*)&(it + 0)->remainders[16 * 0]);
			const uint256_t rems_1 = _mm256_loadu2_m128i((uint128_t*)&(it + 1)->remainders[16 * 1],
														 (uint128_t*)&(it + 0)->remainders[16 * 1]);
			const uint256_t rems_2 = _mm256_loadu2_m128i((uint128_t*)&(it + 1)->remainders[16 * 2],
														 (uint128_t*)&(it + 0)->remainders[16 * 2]);
			const uint256_t rems_3 = _mm256_loadu2_m128i((uint128_t*)&(it + 1)->remainders[16 * 3],
														 (uint128_t*)&(it + 0)->remainders[16 * 3]);
			_mm256_storeu_si256((uint256_t*)(permuted_dt->data() + (32 * 0)), rems_0);
			_mm256_storeu_si256((uint256_t*)(permuted_dt->data() + (32 * 1)), rems_1);
			_mm256_storeu_si256((uint256_t*)(permuted_dt->data() + (32 * 2)), rems_2);
			_mm256_storeu_si256((uint256_t*)(permuted_dt->data() + (32 * 3)), rems_3);
		}
	}



	void full_div_tests::update_div_test_order()
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

		// update the permuted list
		permute_div_tests();
	}

	void full_div_tests::print_div_tests()
	{
	#if analyze_div_tests
		using namespace div_test;

		std::stringstream ss;

		ss << '\n';
		ss << "Prime factor lookup size: " << div_test::detail::prime_factor_lookup.size() << '\n';
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

#if analyze_div_tests
	void full_div_tests::run_div_test_analysis(const uint64_t number)
	{
		using namespace div_test;

		const auto div_test_pred = [](const auto& a, const auto& b)
		{
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
	}
#else
	void full_div_tests::run_div_test_analysis(__attribute__((unused)) const uint64_t number) {}
#endif
}
