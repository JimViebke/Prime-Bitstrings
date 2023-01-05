
#include "multibase_div_tests.hpp"
#include "../util/simd.hpp"

#if analyze_div_tests
#include <sstream>
#endif

namespace mbp::div_test
{
	namespace detail
	{
		consteval std::vector<div_test_t> generate_div_tests_impl()
		{
			std::vector<uncompressed_div_test_t> uncompressed_dts;

			for (size_t i = 1; i < n_of_primes; ++i) // for each small prime starting from 3
			{
				for (size_t base = 3; base <= up_to_base; ++base) // for each base 3..n
				{
					const auto p = small_primes_lookup[i];

					// Skip always-composite cases
					if (base % p == 0) continue;

					// Always suppress hardcoded div tests

					// Implemented with bit pattern filter
					if (base == 3 && p == 5) continue;

					// Hardcoded div tests with 4 remainders
					if (base == 5 && p == 13) continue;
					if (base == 8 && p == 13) continue;
					if (base == 4 && p == 17) continue;

					// 4 remainders, standalone
					if (base == 13 && p == 17) continue;

					// 3 remainders
					if (base == 4 && p == 7) continue;
					if (base == 3 && p == 13) continue;
					if (base == 9 && p == 13) continue;

					// second set with 4 remainders
					if (base == 12 && p == 29) continue;
					if (base == 6 && p == 37) continue;

					// 6 remainders
					if (base == 3 && p == 7) continue;
					if (base == 5 && p == 7) continue;
					if (base == 4 && p == 13) continue;
					if (base == 10 && p == 13) continue;

					// 5 remainders
					if (base == 3 && p == 11) continue;
					if (base == 4 && p == 11) continue;
					if (base == 5 && p == 11) continue;
					if (base == 9 && p == 11) continue;

					// 16 remainders
					if (base == 3 && p == 17) continue;
					if (base == 5 && p == 17) continue;
					if (base == 6 && p == 17) continue;
					if (base == 7 && p == 17) continue;
					if (base == 10 && p == 17) continue;
					if (base == 11 && p == 17) continue;
					if (base == 12 && p == 17) continue;

					// 10 remainders
					if (base == 6 && p == 11) continue;
					if (base == 7 && p == 11) continue;
					if (base == 8 && p == 11) continue;

					// 12 remainders
					if (base == 6 && p == 13) continue;
					if (base == 7 && p == 13) continue;
					if (base == 11 && p == 13) continue;

					// 8 remainders
					if (base == 8 && p == 17) continue;
					if (base == 9 && p == 17) continue;

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

			// Div tests are periodically re-ordered during runtime based on performance.
			// To reflect real-world performance in benchmarks, pre-order the div tests based
			// on cached performance data (the hitcount of each test).
			// Otherwise, start with the default ordering of primes, low to high.
			if constexpr (benchmark_mode)
			{
				for (auto& dt : uncompressed_dts)
				{
					dt.hits = cached_hitcount_for(dt.base, small_primes_lookup[dt.prime_idx]);
				}

				std::sort(uncompressed_dts.begin(), uncompressed_dts.end(),
						  [](const auto& a, const auto& b) { return a.hits > b.hits; });
			}

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
	}

	constexpr size_t div_tests_size = detail::generate_div_tests_impl().size();
	consteval std::array<div_test_t, div_tests_size> generate_div_tests()
	{
		std::array<div_test_t, div_tests_size> div_tests{};
		const auto x = detail::generate_div_tests_impl();
		std::copy(x.begin(), x.end(), div_tests.begin());
		return div_tests;
	}

	using div_tests_t = std::array<div_test::div_test_t, div_test::div_tests_size>;
	div_tests_t div_tests = generate_div_tests(); // intellisense false positive

	namespace detail
	{
		consteval size_t calculate_prime_factor_lookup_size()
		{
			constexpr div_tests_t div_tests_constexpr = generate_div_tests();

			size_t largest_sum = 0;

			// Calculate the largest possible sum of remainders of a number with every bit set
			for (const auto& div_test : div_tests_constexpr)
			{
				const size_t sum = std::accumulate(div_test.remainders.begin(), div_test.remainders.end(), size_t(0));

				if (sum > largest_sum)
					largest_sum = sum;
			}

			// + 1 so "lookup[sum]" is always in range
			return largest_sum + 1;
		}

		constexpr size_t prime_factor_lookup_size = calculate_prime_factor_lookup_size(); // intellisense false positive

		// Faster version
		std::vector<prime_lookup_t> build_prime_factor_lookup_old()
		{
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
	}

	namespace detail
	{
		std::array<std::vector<uint8_t>, n_of_primes> build_indivisible_lookup()
		{
			std::array<std::vector<uint8_t>, n_of_primes> lookup;

			// start from 5, at idx 2
			for (size_t prime_idx = 2; prime_idx < n_of_primes; ++prime_idx)
			{
				// calculate the largest possible sum of remainders, for any
				// div test using this prime

				size_t largest_sum = 0;

				for (const auto& div_test : div_tests)
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

					for (size_t base = 3; base < up_to_base; ++base)
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
		std::array<std::vector<uint8_t>, n_of_primes> indivisible_by = build_indivisible_lookup();
	}



	template<bool on_fast_path>
	size_t* branchless_div_tests(size_t* const candidates_begin,
								 size_t* const candidates_end,
								 const size_t n_of_tests)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		size_t* shrinking_end = candidates_end;

		static constexpr uint256_t static_shuffle_mask_byte_0{ .m256i_u64{
			0x0000000000000000, 0x0808080808080808, 0x0000000000000000, 0x0808080808080808 } };
		static constexpr uint256_t static_shuffle_mask_byte_1{ .m256i_u64{
			0x0101010101010101, 0x0909090909090909, 0x0101010101010101, 0x0909090909090909 } };
		static constexpr uint256_t static_shuffle_mask_byte_2{ .m256i_u64{
			0x0202020202020202, 0x0A0A0A0A0A0A0A0A, 0x0202020202020202, 0x0A0A0A0A0A0A0A0A } };
		static constexpr uint256_t static_shuffle_mask_byte_3{ .m256i_u64{
			0x0303030303030303, 0x0B0B0B0B0B0B0B0B, 0x0303030303030303, 0x0B0B0B0B0B0B0B0B } };

		const uint256_t shuffle_mask_byte_0 = _mm256_load_si256(&static_shuffle_mask_byte_0);
		const uint256_t shuffle_mask_byte_1 = _mm256_load_si256(&static_shuffle_mask_byte_1);
		const uint256_t shuffle_mask_byte_2 = _mm256_load_si256(&static_shuffle_mask_byte_2);
		const uint256_t shuffle_mask_byte_3 = _mm256_load_si256(&static_shuffle_mask_byte_3);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		for (size_t i = 0; i < n_of_tests; ++i)
		{
			div_test_t& div_test = div_tests[i];

			const size_t* input = candidates_begin;
			size_t* output = candidates_begin;

			size_t upper_bits = *input & upper_bits_mask;
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
				const size_t* const rounded_end = input + (n_of_candidates - (n_of_candidates % 4));

				// load four candidates
				uint256_t candidates = _mm256_loadu_si256((uint256_t*)input);

				alignas(32) uint64_t sums[4]{};

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

			size_t number = *input; // load one iteration ahead

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

	template size_t* branchless_div_tests<true>(size_t* const, size_t* const, const size_t);
	template size_t* branchless_div_tests<false>(size_t* const, size_t* const, const size_t);

	template<bool on_fast_path>
	size_t* branching_div_tests(size_t* input,
								const size_t* const candidates_end,
								const size_t start_offset)
	{
		using namespace div_test;

		static_assert(sizeof(remainder_t) == 1);

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		size_t upper_bits = *input & upper_bits_mask;
		uint256_t mask_upper = util::expand_bits_to_bytes(*input >> 32);

		size_t* output = input;

		div_test::div_test_t* div_tests_start = div_tests.data() + start_offset;
		const auto* const div_tests_end = div_tests.data() + div_tests.size();

		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;

			// always write
			*output = number;

			// Convert each half of the number to a 32-byte bitmask

			if constexpr (!on_fast_path)
			{
				// Upper 32 bits only change every 4B integers - re-calculate when stale
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;
					mask_upper = util::expand_bits_to_bytes(number >> 32);
				}
			}

			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));

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
				ymm0 = _mm256_loadu_si256((uint256_t*)(&div_test.remainders[0] + sizeof(div_test_t)));
				ymm1 = _mm256_loadu_si256((uint256_t*)(&div_test.remainders[32] + sizeof(div_test_t)));

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

	template size_t* branching_div_tests<true>(size_t*, const size_t* const, const size_t);
	template size_t* branching_div_tests<false>(size_t*, const size_t* const, const size_t);



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
}
