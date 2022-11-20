#pragma once

#include "config.hpp"
#include "math/math.hpp"
#include "trial_division/types.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/simd.hpp"

namespace mbp
{
	template<bool on_fast_path>
	tests_are_inlined size_t* branchless_div_tests(size_t* const candidates_begin,
												   size_t* const candidates_end,
												   div_test::div_test_t* div_tests,
												   const size_t n_of_tests)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		size_t* shrinking_end = candidates_end;

		for (size_t i = 0; i < n_of_tests; ++i)
		{
			div_test_t& div_test = *(div_tests + i);

			const uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)&div_test.remainders[0]);
			const uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)&div_test.remainders[32]);

			const size_t* input = candidates_begin;
			size_t* output = candidates_begin;

			size_t upper_bits = *input & upper_bits_mask;
			const prime_lookup_t* pf_lookup_preinc = prime_factor_lookup.data() +
				util::vcl_hadd_x(_mm256_and_si256(util::expand_bits_to_bytes(upper_bits >> 32), ymm1));

			if constexpr (on_fast_path)
			{
				const uint256_t shuffle_mask_lo = _mm256_set_epi64x(0x0909090909090909, 0x0808080808080808, 0x0101010101010101, 0x0000000000000000);
				const uint256_t shuffle_mask_hi = _mm256_set_epi64x(0x0B0B0B0B0B0B0B0B, 0x0A0A0A0A0A0A0A0A, 0x0303030303030303, 0x0202020202020202);
				const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

				// We're testing candidates against the same div test.
				// Load the lowest 16 bytes of a div test into high and low lanes of a register
				// Load the next 16 bytes of a div test into high and low lanes of a register
				const uint128_t rems_0_15 = _mm_loadu_si128((uint128_t*)&div_test.remainders[0]);
				const uint256_t rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(rems_0_15), rems_0_15, 1);
				const uint128_t rems_16_31 = _mm_loadu_si128((uint128_t*)&div_test.remainders[16]);
				const uint256_t rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(rems_16_31), rems_16_31, 1);

				// calculate the number of candidates, rounded down to the nearest 2
				const size_t n_of_candidates = shrinking_end - input;
				const size_t* const rounded_end = input + (n_of_candidates - (n_of_candidates % 2));

				// load two candidates
				uint128_t xmm_candidates = _mm_loadu_si128((uint128_t*)input);
				uint256_t candidates = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidates), xmm_candidates, 1);

				for (; input < rounded_end; )
				{
					// convert bits to bytes
					uint256_t candidates_lo = _mm256_shuffle_epi8(candidates, shuffle_mask_lo);
					uint256_t candidates_hi = _mm256_shuffle_epi8(candidates, shuffle_mask_hi);

					// load ahead
					xmm_candidates = _mm_loadu_si128((uint128_t*)(input + 2));
					candidates = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidates), xmm_candidates, 1);

					// convert bits to bytes, cont'd
					candidates_lo = _mm256_andnot_si256(candidates_lo, and_mask);
					candidates_hi = _mm256_andnot_si256(candidates_hi, and_mask);
					candidates_lo = _mm256_cmpeq_epi8(candidates_lo, _mm256_setzero_si256());
					candidates_hi = _mm256_cmpeq_epi8(candidates_hi, _mm256_setzero_si256());

					// mask to select remainders
					candidates_lo = _mm256_and_si256(candidates_lo, rems_lo);
					candidates_hi = _mm256_and_si256(candidates_hi, rems_hi);

					// vertically sum low and high remainders (bytes [0, 1, 0, 1] + [2, 3, 2, 3])
					uint256_t remainders = _mm256_add_epi8(candidates_lo, candidates_hi);
					// h-sum 4 sets of 8 consecutive bytes
					remainders = _mm256_sad_epu8(remainders, _mm256_setzero_si256());
					// final h-sum
					remainders = _mm256_add_epi16(remainders, _mm256_shuffle_epi32(remainders, 2));

					const auto lookup_a = pf_lookup_preinc[remainders.m256i_u32[0]];
					const auto lookup_b = pf_lookup_preinc[remainders.m256i_u32[4]];

					// load and increment
					const size_t candidate_a = *input;
					const size_t candidate_b = *(input + 1);
					input += 2;

					*output = candidate_a; // always write
					// branchless conditional increment
					if ((lookup_a & (1 << div_test.prime_idx)) == 0) ++output;

					*output = candidate_b;
					if ((lookup_b & (1 << div_test.prime_idx)) == 0) ++output;
				}

				// if there is a remaining (odd) candidate, the loop below runs once and handles it

			} // end if (on_fast_path)

			size_t number = *input; // load one iteration ahead

			for (; input < shrinking_end; )
			{
				if constexpr (!on_fast_path)
				{
					// Upper 32 bits only change every 4B integers - re-calculate when stale
					if ((number & upper_bits_mask) != upper_bits)
					{
						upper_bits = number & upper_bits_mask;
						pf_lookup_preinc = prime_factor_lookup.data() +
							util::vcl_hadd_x(_mm256_and_si256(util::expand_bits_to_bytes(upper_bits >> 32), ymm1));
					}
				}

				const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));

				*output = number; // always write
				number = *++input; // load one iteration ahead

				const uint256_t rems_lower = _mm256_and_si256(mask_lower, ymm0);
				const auto rem = util::vcl_hadd_x(rems_lower);

				// increment only if we haven't found a prime factor yet
				size_t lookup = pf_lookup_preinc[rem] >> div_test.prime_idx;
				output += ~lookup & 0b1;
			}

			div_test.hits += uint32_t(input - output);

			// output is one past the last element; use it as the new end
			shrinking_end = output;
		}

		return shrinking_end;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* multibase_div_tests(size_t* input,
												  const size_t* const candidates_end,
												  div_test::div_test_t* div_tests_start)
	{
		using namespace div_test;

		static_assert(sizeof(remainder_t) == 1);

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		size_t upper_bits = *input & upper_bits_mask;
		uint256_t mask_upper = util::expand_bits_to_bytes(*input >> 32);

		size_t* output = input;
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

}
