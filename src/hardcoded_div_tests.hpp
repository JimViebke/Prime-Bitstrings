#pragma once

#include "config.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/simd.hpp"

namespace mbp
{
	namespace detail
	{
		template<size_t base, size_t place_value_start, size_t prime>
		consteval uint256_t build_shuffle_lookup_impl()
		{
			uint256_t lookup{ .m256i_u8 = { 0 } };
			for (size_t i = 0; i < 16; ++i)
			{
				if (i & 0b0001) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 0, prime>::rem;
				if (i & 0b0010) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 1, prime>::rem;
				if (i & 0b0100) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 2, prime>::rem;
				if (i & 0b1000) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 3, prime>::rem;
				lookup.m256i_u8[i + 16] = lookup.m256i_u8[i]; // duplicate into upper lane
			}
			return lookup;
		}

		template<size_t base, size_t prime>
		consteval uint256_t build_3rem_shuffle_lookup()
		{
			uint256_t lookup{ .m256i_u8 = { 0 } };
			for (size_t i = 0; i < 8; ++i)
			{
				if (i & 0b001) lookup.m256i_u8[i] += pow_mod<base, 0, prime>::rem;
				if (i & 0b010) lookup.m256i_u8[i] += pow_mod<base, 1, prime>::rem;
				if (i & 0b100) lookup.m256i_u8[i] += pow_mod<base, 2, prime>::rem;
				lookup.m256i_u8[i + 16] = lookup.m256i_u8[i]; // duplicate into upper lane
			}
			return lookup;
		}

		template<size_t base, size_t prime>
		consteval uint256_t build_4rem_shuffle_lookup()
		{
			return build_shuffle_lookup_impl<base, 0, prime>();
		}

		template<size_t base, size_t prime>
		consteval uint256_t build_5rem_shuffle_lookup()
		{
			return build_shuffle_lookup_impl<base, 1, prime>();
		}

		template<size_t base, size_t prime>
		consteval uint256_t build_8rem_shuffle_lookup_lo_nybble()
		{
			return build_shuffle_lookup_impl<base, 0, prime>();
		}

		template<size_t base, size_t prime>
		consteval uint256_t build_8rem_shuffle_lookup_hi_nybble()
		{
			return build_shuffle_lookup_impl<base, 4, prime>();
		}
	}

	template<bool on_fast_path>
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

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) & 0b111);

			constexpr static uint256_t static_nybble_lookup_a = mbp::detail::build_4rem_shuffle_lookup<3, 5>();
			constexpr static uint256_t static_nybble_lookup_b = mbp::detail::build_4rem_shuffle_lookup<5, 13>();
			constexpr static uint256_t static_nybble_lookup_c = mbp::detail::build_4rem_shuffle_lookup<8, 13>();
			constexpr static uint256_t static_nybble_lookup_d = mbp::detail::build_4rem_shuffle_lookup<4, 17>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup_a = _mm256_loadu_si256(&static_nybble_lookup_a);
			const uint256_t nybble_lookup_b = _mm256_loadu_si256(&static_nybble_lookup_b);
			const uint256_t nybble_lookup_c = _mm256_loadu_si256(&static_nybble_lookup_c);
			const uint256_t nybble_lookup_d = _mm256_loadu_si256(&static_nybble_lookup_d);

			// calculate upper sum of remainders
			const size_t number_upper = (*input) & (size_t(-1) << 32);
			const size_t upc_0 = pop_count(number_upper & (bitmask << 0));
			size_t b3m5_urem = upc_0;
			size_t b8m13_urem = upc_0;
			size_t b5m13_urem = upc_0;
			size_t b4m17_urem = upc_0;
			const size_t upc_1 = pop_count(number_upper & (bitmask << 1));
			b3m5_urem += upc_1 * pow_mod<3, 1, 5>::rem;
			b5m13_urem += upc_1 * pow_mod<5, 1, 13>::rem;
			b8m13_urem += upc_1 * pow_mod<8, 1, 13>::rem;
			b4m17_urem += upc_1 * pow_mod<4, 1, 17>::rem;
			const size_t upc_2 = pop_count(number_upper & (bitmask << 2));
			b3m5_urem += upc_2 * pow_mod<3, 2, 5>::rem;
			b5m13_urem += upc_2 * pow_mod<5, 2, 13>::rem;
			b8m13_urem += upc_2 * pow_mod<8, 2, 13>::rem;
			b4m17_urem += upc_2 * pow_mod<4, 2, 17>::rem;
			const size_t upc_3 = pop_count(number_upper & (bitmask << 3));
			b3m5_urem += upc_3 * pow_mod<3, 3, 5>::rem;
			b5m13_urem += upc_3 * pow_mod<5, 3, 13>::rem;
			b8m13_urem += upc_3 * pow_mod<8, 3, 13>::rem;
			b4m17_urem += upc_3 * pow_mod<4, 3, 17>::rem;
			const uint8_t* const indivisible_by_5_ptr = indivisible_by[get_prime_index<5>::idx].data() + b3m5_urem;
			const uint8_t* const indivisible_by_13_ptr_b5 = indivisible_by[get_prime_index<13>::idx].data() + b5m13_urem;
			const uint8_t* const indivisible_by_13_ptr_b8 = indivisible_by[get_prime_index<13>::idx].data() + b8m13_urem;
			const uint8_t* const indivisible_by_17_ptr = indivisible_by[get_prime_index<17>::idx].data() + b4m17_urem;

			alignas(32) uint16_t sums_ab[16]{};
			alignas(32) uint16_t sums_cd[16]{};

			// run vector instructions one iteration ahead
			{
				// load 8 candidates into two registers
				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));

				// we only want the bottom 32 bits of our 8 candidates
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				// select the high and low halves of each byte
				const uint256_t nybbles_lo = _mm256_and_si256(ymm0, nybble_mask);
				const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);
				const uint256_t candidates_0426 = _mm256_unpacklo_epi32(nybbles_lo, nybbles_hi);
				const uint256_t candidates_1537 = _mm256_unpackhi_epi32(nybbles_lo, nybbles_hi);

				// replace the bits of each nybble with their remainder, then sum remainders
				uint256_t rems_0426 = _mm256_shuffle_epi8(nybble_lookup_a, candidates_0426);
				uint256_t rems_1537 = _mm256_shuffle_epi8(nybble_lookup_a, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_a = _mm256_packus_epi32(rems_0426, rems_1537);

				rems_0426 = _mm256_shuffle_epi8(nybble_lookup_b, candidates_0426);
				rems_1537 = _mm256_shuffle_epi8(nybble_lookup_b, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_b = _mm256_packus_epi32(rems_0426, rems_1537);

				const uint256_t sums_0415a_0415b_2637a_2637b = _mm256_packus_epi32(sums_04152637_a, sums_04152637_b);

				rems_0426 = _mm256_shuffle_epi8(nybble_lookup_c, candidates_0426);
				rems_1537 = _mm256_shuffle_epi8(nybble_lookup_c, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_c = _mm256_packus_epi32(rems_0426, rems_1537);

				rems_0426 = _mm256_shuffle_epi8(nybble_lookup_d, candidates_0426);
				rems_1537 = _mm256_shuffle_epi8(nybble_lookup_d, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_d = _mm256_packus_epi32(rems_0426, rems_1537);

				const uint256_t sums_0415c_0415d_2637c_2637d = _mm256_packus_epi32(sums_04152637_c, sums_04152637_d);

				// store results on the stack
				_mm256_store_si256((uint256_t*)sums_ab, sums_0415a_0415b_2637a_2637b);
				_mm256_store_si256((uint256_t*)sums_cd, sums_0415c_0415d_2637c_2637d);
			}

			// load two iterations ahead
			uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 8));
			uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 12));

			for (; input < candidates_end_rounded; )
			{
				// run vector instructions one iteration ahead

				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				const uint256_t nybbles_lo = _mm256_and_si256(ymm0, nybble_mask);
				const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				const uint256_t candidates_0426 = _mm256_unpacklo_epi32(nybbles_lo, nybbles_hi);
				const uint256_t candidates_1537 = _mm256_unpackhi_epi32(nybbles_lo, nybbles_hi);

				uint256_t rems_0426 = _mm256_shuffle_epi8(nybble_lookup_a, candidates_0426);
				uint256_t rems_1537 = _mm256_shuffle_epi8(nybble_lookup_a, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_a = _mm256_packus_epi32(rems_0426, rems_1537);

				rems_0426 = _mm256_shuffle_epi8(nybble_lookup_b, candidates_0426);
				rems_1537 = _mm256_shuffle_epi8(nybble_lookup_b, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_b = _mm256_packus_epi32(rems_0426, rems_1537);

				const uint256_t sums_0415a_0415b_2637a_2637b = _mm256_packus_epi32(sums_04152637_a, sums_04152637_b);

				rems_0426 = _mm256_shuffle_epi8(nybble_lookup_c, candidates_0426);
				rems_1537 = _mm256_shuffle_epi8(nybble_lookup_c, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_c = _mm256_packus_epi32(rems_0426, rems_1537);

				rems_0426 = _mm256_shuffle_epi8(nybble_lookup_d, candidates_0426);
				rems_1537 = _mm256_shuffle_epi8(nybble_lookup_d, candidates_1537);
				rems_0426 = _mm256_sad_epu8(rems_0426, _mm256_setzero_si256());
				rems_1537 = _mm256_sad_epu8(rems_1537, _mm256_setzero_si256());
				const uint256_t sums_04152637_d = _mm256_packus_epi32(rems_0426, rems_1537);

				const uint256_t sums_0415c_0415d_2637c_2637d = _mm256_packus_epi32(sums_04152637_c, sums_04152637_d);

				// load two iterations ahead
				ymm0 = _mm256_loadu_si256((uint256_t*)(input + 16));
				ymm1 = _mm256_loadu_si256((uint256_t*)(input + 20));

				// Only advance the pointer if the number is still a candidate, that is,
				// it is not known to be divisible for the above tests

				const size_t inc_0 = indivisible_by_5_ptr[sums_ab[0]]
					& indivisible_by_13_ptr_b5[sums_ab[0 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[0]]
					& indivisible_by_17_ptr[sums_cd[0 + 4]];
				*output = *input++; // always copy
				output = (uint64_t*)(((uint8_t*)output) + inc_0); // branchless conditional increment

				const size_t inc_1 = indivisible_by_5_ptr[sums_ab[2]]
					& indivisible_by_13_ptr_b5[sums_ab[2 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[2]]
					& indivisible_by_17_ptr[sums_cd[2 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_1);

				const size_t inc_2 = indivisible_by_5_ptr[sums_ab[8]]
					& indivisible_by_13_ptr_b5[sums_ab[8 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[8]]
					& indivisible_by_17_ptr[sums_cd[8 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_2);

				const size_t inc_3 = indivisible_by_5_ptr[sums_ab[10]]
					& indivisible_by_13_ptr_b5[sums_ab[10 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[10]]
					& indivisible_by_17_ptr[sums_cd[10 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_3);

				const size_t inc_4 = indivisible_by_5_ptr[sums_ab[1]]
					& indivisible_by_13_ptr_b5[sums_ab[1 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[1]]
					& indivisible_by_17_ptr[sums_cd[1 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_4);

				const size_t inc_5 = indivisible_by_5_ptr[sums_ab[3]]
					& indivisible_by_13_ptr_b5[sums_ab[3 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[3]]
					& indivisible_by_17_ptr[sums_cd[3 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_5);

				const size_t inc_6 = indivisible_by_5_ptr[sums_ab[9]]
					& indivisible_by_13_ptr_b5[sums_ab[9 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[9]]
					& indivisible_by_17_ptr[sums_cd[9 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_6);

				const size_t inc_7 = indivisible_by_5_ptr[sums_ab[11]]
					& indivisible_by_13_ptr_b5[sums_ab[11 + 4]]
					& indivisible_by_13_ptr_b8[sums_cd[11]]
					& indivisible_by_17_ptr[sums_cd[11 + 4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_7);

				// store on the stack for the next iteration
				_mm256_store_si256((uint256_t*)sums_ab, sums_0415a_0415b_2637a_2637b);
				_mm256_store_si256((uint256_t*)sums_cd, sums_0415c_0415d_2637c_2637d);
			}
		}

		const uint8_t* const indivisible_by_5 = indivisible_by[get_prime_index<5>::idx].data();
		const uint8_t* const indivisible_by_13 = indivisible_by[get_prime_index<13>::idx].data();
		const uint8_t* const indivisible_by_17 = indivisible_by[get_prime_index<17>::idx].data();

		// handle any remaining elements
		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;
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

			const size_t inc = indivisible_by_5[b3m5_rem]
				& indivisible_by_13[b5m13_rem]
				& indivisible_by_13[b8m13_rem]
				& indivisible_by_17[b4m17_rem];

			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
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

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static uint256_t static_nybble_lookup_b4m7 = mbp::detail::build_3rem_shuffle_lookup<4, 7>();
			constexpr static uint256_t static_nybble_lookup_b3m13 = mbp::detail::build_3rem_shuffle_lookup<3, 13>();
			constexpr static uint256_t static_nybble_lookup_b9m13 = mbp::detail::build_3rem_shuffle_lookup<9, 13>();

			constexpr uint64_t pdep_mask_lo = 0x0707070707070707; // first 24 bits
			constexpr uint64_t pdep_mask_hi = 0x0000000000030707; // next 8 (24 + 8 == 32)

			constexpr size_t upper_bits_mask = size_t(-1) << 32;
			uint64_t upper_bits = (*input) & upper_bits_mask;
			const size_t pc_0 = pop_count(upper_bits & (bitmask << 0));
			size_t upper_sum_b4m7 = pc_0;
			size_t upper_sum_b3m13 = pc_0;
			size_t upper_sum_b9m13 = pc_0;
			const size_t pc_1 = pop_count(upper_bits & (bitmask << 1));
			upper_sum_b4m7 += pc_1 * pow_mod<4, 1, 7>::rem;
			upper_sum_b3m13 += pc_1 * pow_mod<3, 1, 13>::rem;
			upper_sum_b9m13 += pc_1 * pow_mod<9, 1, 13>::rem;
			const size_t pc_2 = pop_count(upper_bits & (bitmask << 2));
			upper_sum_b4m7 += pc_2 * pow_mod<4, 2, 7>::rem;
			upper_sum_b3m13 += pc_2 * pow_mod<3, 2, 13>::rem;
			upper_sum_b9m13 += pc_2 * pow_mod<9, 2, 13>::rem;
			const uint8_t* const indivisible_b4m7 = indivisible_by[get_prime_index<7>::idx].data() + upper_sum_b4m7;
			const uint8_t* const indivisible_b3m13 = indivisible_by[get_prime_index<13>::idx].data() + upper_sum_b3m13;
			const uint8_t* const indivisible_b9m13 = indivisible_by[get_prime_index<13>::idx].data() + upper_sum_b9m13;

			const uint256_t nybble_lookup_b4m7 = _mm256_loadu_si256(&static_nybble_lookup_b4m7);
			const uint256_t nybble_lookup_b3m13 = _mm256_loadu_si256(&static_nybble_lookup_b3m13);
			const uint256_t nybble_lookup_b9m13 = _mm256_loadu_si256(&static_nybble_lookup_b9m13);

			const uint64_t* const rounded_end = candidates_end - ((candidates_end - input) & 0b11);

			alignas(32) uint64_t pdep_candidates[8]{};
			alignas(32) uint16_t rems[16]{};

			// run pdep and vector instructions one iteration ahead
			{
				// move each set of three bits to its own byte
				uint64_t number = *(input + 0);
				pdep_candidates[0] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[4] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 1);
				pdep_candidates[1] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[5] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 2);
				pdep_candidates[2] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[6] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 3);
				pdep_candidates[3] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[7] = _pdep_u64(number >> 24, pdep_mask_hi);

				const uint256_t candidates_lo = _mm256_loadu_si256((uint256_t*)(pdep_candidates + 0));
				const uint256_t candidates_hi = _mm256_loadu_si256((uint256_t*)(pdep_candidates + 4));

				// replace each set of three bits with its sum of remainders
				const uint256_t rems_a_lo = _mm256_shuffle_epi8(nybble_lookup_b4m7, candidates_lo);
				const uint256_t rems_a_hi = _mm256_shuffle_epi8(nybble_lookup_b4m7, candidates_hi);
				const uint256_t rems_b_lo = _mm256_shuffle_epi8(nybble_lookup_b3m13, candidates_lo);
				const uint256_t rems_b_hi = _mm256_shuffle_epi8(nybble_lookup_b3m13, candidates_hi);
				const uint256_t rems_c_lo = _mm256_shuffle_epi8(nybble_lookup_b9m13, candidates_lo);
				const uint256_t rems_c_hi = _mm256_shuffle_epi8(nybble_lookup_b9m13, candidates_hi);

				// vertically sum remainders
				uint256_t rems_a = _mm256_add_epi8(rems_a_lo, rems_a_hi);
				uint256_t rems_b = _mm256_add_epi8(rems_b_lo, rems_b_hi);
				uint256_t rems_c = _mm256_add_epi8(rems_c_lo, rems_c_hi);
				// h-sum remainders
				rems_a = _mm256_sad_epu8(rems_a, _mm256_setzero_si256());
				rems_b = _mm256_sad_epu8(rems_b, _mm256_setzero_si256());
				rems_c = _mm256_sad_epu8(rems_c, _mm256_setzero_si256());

				uint256_t rems_ab = _mm256_packus_epi32(rems_a, rems_b);
				uint256_t rems_cx = _mm256_packus_epi32(rems_c, _mm256_setzero_si256());
				uint256_t rems_abcx = _mm256_packus_epi32(rems_ab, rems_cx);

				// store on the stack
				_mm256_storeu_si256((uint256_t*)rems, rems_abcx);
			}

			// run pdep instructions two iterations ahead
			{
				uint64_t number = *(input + 4 + 0);
				pdep_candidates[0] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[4] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 4 + 1);
				pdep_candidates[1] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[5] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 4 + 2);
				pdep_candidates[2] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[6] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 4 + 3);
				pdep_candidates[3] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[7] = _pdep_u64(number >> 24, pdep_mask_hi);
			}

			uint256_t candidates_lo = _mm256_loadu_si256((uint256_t*)(pdep_candidates + 0));
			uint256_t candidates_hi = _mm256_loadu_si256((uint256_t*)(pdep_candidates + 4));

			for (; input < rounded_end; )
			{
				// load four candidates two iterations ahead
				uint64_t number = *(input + 8 + 0);
				pdep_candidates[0] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[4] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 8 + 1);
				pdep_candidates[1] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[5] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 8 + 2);
				pdep_candidates[2] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[6] = _pdep_u64(number >> 24, pdep_mask_hi);
				number = *(input + 8 + 3);
				pdep_candidates[3] = _pdep_u64(number, pdep_mask_lo);
				pdep_candidates[7] = _pdep_u64(number >> 24, pdep_mask_hi);

				// run vector instructions one iteration ahead
				const uint256_t rems_a_lo = _mm256_shuffle_epi8(nybble_lookup_b4m7, candidates_lo);
				const uint256_t rems_a_hi = _mm256_shuffle_epi8(nybble_lookup_b4m7, candidates_hi);
				const uint256_t rems_b_lo = _mm256_shuffle_epi8(nybble_lookup_b3m13, candidates_lo);
				const uint256_t rems_b_hi = _mm256_shuffle_epi8(nybble_lookup_b3m13, candidates_hi);
				const uint256_t rems_c_lo = _mm256_shuffle_epi8(nybble_lookup_b9m13, candidates_lo);
				const uint256_t rems_c_hi = _mm256_shuffle_epi8(nybble_lookup_b9m13, candidates_hi);

				// load for the next iteration
				candidates_lo = _mm256_loadu_si256((uint256_t*)(pdep_candidates + 0));
				candidates_hi = _mm256_loadu_si256((uint256_t*)(pdep_candidates + 4));

				uint256_t rems_a = _mm256_add_epi8(rems_a_lo, rems_a_hi);
				uint256_t rems_b = _mm256_add_epi8(rems_b_lo, rems_b_hi);
				uint256_t rems_c = _mm256_add_epi8(rems_c_lo, rems_c_hi);
				rems_a = _mm256_sad_epu8(rems_a, _mm256_setzero_si256());
				rems_b = _mm256_sad_epu8(rems_b, _mm256_setzero_si256());
				rems_c = _mm256_sad_epu8(rems_c, _mm256_setzero_si256());

				uint256_t rems_ab = _mm256_packus_epi32(rems_a, rems_b);
				uint256_t rems_cx = _mm256_packus_epi32(rems_c, _mm256_setzero_si256());
				uint256_t rems_abcx = _mm256_packus_epi32(rems_ab, rems_cx);

				const size_t inc_0 = indivisible_b4m7[rems[0 + 0]]
					& indivisible_b3m13[rems[2 + 0]]
					& indivisible_b9m13[rems[4 + 0]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_0);

				const size_t inc_1 = indivisible_b4m7[rems[0 + 1]]
					& indivisible_b3m13[rems[2 + 1]]
					& indivisible_b9m13[rems[4 + 1]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_1);

				const size_t inc_2 = indivisible_b4m7[rems[0 + 8]]
					& indivisible_b3m13[rems[2 + 8]]
					& indivisible_b9m13[rems[4 + 8]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_2);

				const size_t inc_3 = indivisible_b4m7[rems[0 + 9]]
					& indivisible_b3m13[rems[2 + 9]]
					& indivisible_b9m13[rems[4 + 9]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_3);

				// store above results on the stack for the next iteration
				_mm256_storeu_si256((uint256_t*)rems, rems_abcx);
			}
		}

		const uint8_t* const indivisible_by_7 = indivisible_by[get_prime_index<7>::idx].data();
		const uint8_t* const indivisible_by_13 = indivisible_by[get_prime_index<13>::idx].data();

		size_t number = *input;

		for (; input < candidates_end; )
		{
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

			*output = number; // always write
			number = *++input; // load ahead

			// Only advance the pointer if the number is still a candidate
			size_t inc = indivisible_by_7[b4_m7_rem]
				& indivisible_by_13[b3_m13_rem]
				& indivisible_by_13[b9_m13_rem];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
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

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		constexpr static uint256_t static_rems_byte_0 = uint256_t{ .m256i_u8{
			1, 3, 2, 6, 4, 5, 1, 3, 1, 5, 4, 6, 2, 3, 1, 5, 1, 4, 3, 12, 9, 10, 1, 4, 1, 10, 9, 12, 3, 4, 1, 10 } };
		constexpr static uint256_t static_rems_byte_1 = uint256_t{ .m256i_u8{
			2, 6, 4, 5, 1, 3, 2, 6, 4, 6, 2, 3, 1, 5, 4, 6, 3, 12, 9, 10, 1, 4, 3, 12, 9, 12, 3, 4, 1, 10, 9, 12 } };
		constexpr static uint256_t static_rems_byte_2 = uint256_t{ .m256i_u8{
			4, 5, 1, 3, 2, 6, 4, 5, 2, 3, 1, 5, 4, 6, 2, 3, 9, 10, 1, 4, 3, 12, 9, 10, 3, 4, 1, 10, 9, 12, 3, 4 } };
		// byte 3 same as byte 0

		size_t* output = input;
		size_t number = *input; // load one iteration ahead

		size_t upper_bits = number & upper_bits_mask;

		const size_t pc_0 = pop_count(upper_bits & (bitmask << 0));
		size_t upper_sum_0 = pc_0;
		size_t upper_sum_1 = pc_0;
		size_t upper_sum_2 = pc_0;
		size_t upper_sum_3 = pc_0;
		const size_t pc_1 = pop_count(upper_bits & (bitmask << 1));
		upper_sum_0 += pc_1 * pow_mod<3, 1, 7>::rem;
		upper_sum_1 += pc_1 * pow_mod<5, 1, 7>::rem;
		upper_sum_2 += pc_1 * pow_mod<4, 1, 13>::rem;
		upper_sum_3 += pc_1 * pow_mod<10, 1, 13>::rem;
		const size_t pc_2 = pop_count(upper_bits & (bitmask << 2));
		upper_sum_0 += pc_2 * pow_mod<3, 2, 7>::rem;
		upper_sum_1 += pc_2 * pow_mod<5, 2, 7>::rem;
		upper_sum_2 += pc_2 * pow_mod<4, 2, 13>::rem;
		upper_sum_3 += pc_2 * pow_mod<10, 2, 13>::rem;
		const size_t pc_3 = pop_count(upper_bits & (bitmask << 3));
		upper_sum_0 += pc_3 * pow_mod<3, 3, 7>::rem;
		upper_sum_1 += pc_3 * pow_mod<5, 3, 7>::rem;
		upper_sum_2 += pc_3 * pow_mod<4, 3, 13>::rem;
		upper_sum_3 += pc_3 * pow_mod<10, 3, 13>::rem;
		const size_t pc_4 = pop_count(upper_bits & (bitmask << 4));
		upper_sum_0 += pc_4 * pow_mod<3, 4, 7>::rem;
		upper_sum_1 += pc_4 * pow_mod<5, 4, 7>::rem;
		upper_sum_2 += pc_4 * pow_mod<4, 4, 13>::rem;
		upper_sum_3 += pc_4 * pow_mod<10, 4, 13>::rem;
		const size_t pc_5 = pop_count(upper_bits & (bitmask << 5));
		upper_sum_0 += pc_5 * pow_mod<3, 5, 7>::rem;
		upper_sum_1 += pc_5 * pow_mod<5, 5, 7>::rem;
		upper_sum_2 += pc_5 * pow_mod<4, 5, 13>::rem;
		upper_sum_3 += pc_5 * pow_mod<10, 5, 13>::rem;

		const uint8_t* indivisible_b3m7 = indivisible_by[get_prime_index<7>::idx].data() + upper_sum_0;
		const uint8_t* indivisible_b5m7 = indivisible_by[get_prime_index<7>::idx].data() + upper_sum_1;
		const uint8_t* indivisible_b4m13 = indivisible_by[get_prime_index<13>::idx].data() + upper_sum_2;
		const uint8_t* indivisible_b10m13 = indivisible_by[get_prime_index<13>::idx].data() + upper_sum_3;

		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		const uint256_t rems_byte_0 = _mm256_loadu_si256(&static_rems_byte_0);
		const uint256_t rems_byte_1 = _mm256_loadu_si256(&static_rems_byte_1);
		const uint256_t rems_byte_2 = _mm256_loadu_si256(&static_rems_byte_2);
		// remainders repeat every three bytes, so instead of having rems_byte_3, re-use rems_byte_0

		const uint256_t shuffle_mask = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, 0x0101010101010101, 0x0000000000000000);

		alignas(32) uint64_t sums[4]{};

		// run vector instructions one iteration ahead
		{
			// this loads more than we need, but VBROADCASTI128 is the cheapest way to load to both lanes
			const uint128_t xmm_candidate = _mm_loadu_si128((uint128_t*)input);
			uint256_t candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

			// convert bits to bytes
			candidate = _mm256_shuffle_epi8(candidate, shuffle_mask);
			candidate = _mm256_andnot_si256(candidate, and_mask);
			candidate = _mm256_cmpeq_epi8(candidate, _mm256_setzero_si256());
			const uint256_t candidate_byte_0 = _mm256_permute4x64_epi64(candidate, 0);
			const uint256_t candidate_byte_1 = _mm256_permute4x64_epi64(candidate, 0b01'01'01'01);
			const uint256_t candidate_byte_2 = _mm256_permute4x64_epi64(candidate, 0b10'10'10'10);
			const uint256_t candidate_byte_3 = _mm256_permute4x64_epi64(candidate, 0b11'11'11'11);

			// select remainders
			const uint256_t sums_0 = _mm256_and_si256(candidate_byte_0, rems_byte_0);
			const uint256_t sums_1 = _mm256_and_si256(candidate_byte_1, rems_byte_1);
			const uint256_t sums_2 = _mm256_and_si256(candidate_byte_2, rems_byte_2);
			const uint256_t sums_3 = _mm256_and_si256(candidate_byte_3, rems_byte_0); // same as 0

			// vertically sum
			const uint256_t sums_01 = _mm256_add_epi8(sums_0, sums_1);
			const uint256_t sums_23 = _mm256_add_epi8(sums_2, sums_3);
			uint256_t sums_0123 = _mm256_add_epi8(sums_01, sums_23);
			// horizontally sum
			sums_0123 = _mm256_sad_epu8(sums_0123, _mm256_setzero_si256());

			// store on the stack
			_mm256_storeu_si256((uint256_t*)sums, sums_0123);
		}

		// this loads more than we need, but VBROADCASTI128 is the cheapest way to load to both lanes
		uint128_t xmm_candidate = _mm_loadu_si128((uint128_t*)(input + 1));
		uint256_t candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

		for (; input < candidates_end; )
		{
			if constexpr (!on_fast_path)
			{
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;

					// recalculate
					indivisible_b3m7 = indivisible_by[get_prime_index<7>::idx].data() + get_upper_sum_of_rems<7, in_base<3>>(number);
					indivisible_b5m7 = indivisible_by[get_prime_index<7>::idx].data() + get_upper_sum_of_rems<7, in_base<5>>(number);
					indivisible_b4m13 = indivisible_by[get_prime_index<13>::idx].data() + get_upper_sum_of_rems<13, in_base<4>>(number);
					indivisible_b10m13 = indivisible_by[get_prime_index<13>::idx].data() + get_upper_sum_of_rems<13, in_base<10>>(number);
				}
			}

			candidate = _mm256_shuffle_epi8(candidate, shuffle_mask);
			candidate = _mm256_andnot_si256(candidate, and_mask);
			candidate = _mm256_cmpeq_epi8(candidate, _mm256_setzero_si256());
			const uint256_t candidate_byte_0 = _mm256_permute4x64_epi64(candidate, 0);
			const uint256_t candidate_byte_1 = _mm256_permute4x64_epi64(candidate, 0b01'01'01'01);
			const uint256_t candidate_byte_2 = _mm256_permute4x64_epi64(candidate, 0b10'10'10'10);
			const uint256_t candidate_byte_3 = _mm256_permute4x64_epi64(candidate, 0b11'11'11'11);

			// load two iterations ahead
			xmm_candidate = _mm_loadu_si128((uint128_t*)(input + 2));
			candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

			const uint256_t sums_0 = _mm256_and_si256(candidate_byte_0, rems_byte_0);
			const uint256_t sums_1 = _mm256_and_si256(candidate_byte_1, rems_byte_1);
			const uint256_t sums_2 = _mm256_and_si256(candidate_byte_2, rems_byte_2);
			const uint256_t sums_3 = _mm256_and_si256(candidate_byte_3, rems_byte_0);

			const uint256_t sums_01 = _mm256_add_epi8(sums_0, sums_1);
			const uint256_t sums_23 = _mm256_add_epi8(sums_2, sums_3);
			uint256_t sums_0123 = _mm256_add_epi8(sums_01, sums_23);
			sums_0123 = _mm256_sad_epu8(sums_0123, _mm256_setzero_si256());

			*output = number; // always write
			number = *++input; // load one iteration ahead

			// Only advance the pointer if the number is still a candidate
			size_t inc = indivisible_b3m7[sums[0]]
				& indivisible_b5m7[sums[1]]
				& indivisible_b4m13[sums[2]]
				& indivisible_b10m13[sums[3]];
			output = (uint64_t*)(((uint8_t*)output) + inc);

			// store on the stack for the next iteration
			_mm256_storeu_si256((uint256_t*)sums, sums_0123);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* div_tests_with_five_rems(size_t* input,
													   const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<3, 11>::val;
		static_assert(bitmask == bitmask_for<4, 11>::val);
		static_assert(bitmask == bitmask_for<5, 11>::val);
		static_assert(bitmask == bitmask_for<9, 11>::val);
		static_assert(period_of<bitmask>::val == 5);

		size_t* output = input;

		const uint8_t* const indivisible_by_11 = indivisible_by[get_prime_index<11>::idx].data();

		if constexpr (on_fast_path)
		{
			constexpr static uint256_t static_nybble_lookup_b3 = mbp::detail::build_5rem_shuffle_lookup<3, 11>();
			constexpr static uint256_t static_nybble_lookup_b4 = mbp::detail::build_5rem_shuffle_lookup<4, 11>();
			constexpr static uint256_t static_nybble_lookup_b5 = mbp::detail::build_5rem_shuffle_lookup<5, 11>();
			constexpr static uint256_t static_nybble_lookup_b9 = mbp::detail::build_5rem_shuffle_lookup<9, 11>();

			constexpr uint64_t upper_bits_mask = uint64_t(-1) << 32;

			const uint64_t upper_bits = *input & upper_bits_mask;
			const size_t pc_0 = pop_count(upper_bits & (bitmask << 0));
			size_t sum_0 = pc_0;
			size_t sum_1 = pc_0;
			size_t sum_2 = pc_0;
			size_t sum_3 = pc_0;
			const size_t pc_1 = pop_count(upper_bits & (bitmask << 1));
			sum_0 += pc_1 * pow_mod<3, 1, 11>::rem;
			sum_1 += pc_1 * pow_mod<4, 1, 11>::rem;
			sum_2 += pc_1 * pow_mod<5, 1, 11>::rem;
			sum_3 += pc_1 * pow_mod<9, 1, 11>::rem;
			const size_t pc_2 = pop_count(upper_bits & (bitmask << 2));
			sum_0 += pc_2 * pow_mod<3, 2, 11>::rem;
			sum_1 += pc_2 * pow_mod<4, 2, 11>::rem;
			sum_2 += pc_2 * pow_mod<5, 2, 11>::rem;
			sum_3 += pc_2 * pow_mod<9, 2, 11>::rem;
			const size_t pc_3 = pop_count(upper_bits & (bitmask << 3));
			sum_0 += pc_3 * pow_mod<3, 3, 11>::rem;
			sum_1 += pc_3 * pow_mod<4, 3, 11>::rem;
			sum_2 += pc_3 * pow_mod<5, 3, 11>::rem;
			sum_3 += pc_3 * pow_mod<9, 3, 11>::rem;
			const size_t pc_4 = pop_count(upper_bits & (bitmask << 4));
			sum_0 += pc_4 * pow_mod<3, 4, 11>::rem;
			sum_1 += pc_4 * pow_mod<4, 4, 11>::rem;
			sum_2 += pc_4 * pow_mod<5, 4, 11>::rem;
			sum_3 += pc_4 * pow_mod<9, 4, 11>::rem;

			const uint8_t* const indivisible_b3 = indivisible_by_11 + sum_0;
			const uint8_t* const indivisible_b4 = indivisible_by_11 + sum_1;
			const uint8_t* const indivisible_b5 = indivisible_by_11 + sum_2;
			const uint8_t* const indivisible_b9 = indivisible_by_11 + sum_3;

			const uint256_t nybble_mask = _mm256_set1_epi8(0b00001111);

			const uint256_t nybble_lookup_b3 = _mm256_loadu_si256(&static_nybble_lookup_b3);
			const uint256_t nybble_lookup_b4 = _mm256_loadu_si256(&static_nybble_lookup_b4);
			const uint256_t nybble_lookup_b5 = _mm256_loadu_si256(&static_nybble_lookup_b5);
			const uint256_t nybble_lookup_b9 = _mm256_loadu_si256(&static_nybble_lookup_b9);

			// Move bits 0,5,10... to the lowest bit of a high nybble, and the other four bits to the (next) low nybble
			// Split this using andmasks once we have four candidates in a vector register
			constexpr uint64_t pdep_mask = 0b1'00011111'00011111'00011111'00011111'00011111'00011111'00010000;
			static_assert(std::popcount(pdep_mask) == 32);

			const uint64_t* const rounded_end = candidates_end - ((candidates_end - input) & 0b11);

			alignas(32) uint64_t pdep_candidates[4]{};
			alignas(32) uint16_t sums[16]{};

			// run pdep and vector instructions one iteration ahead
			{
				uint64_t number = *(input + 0);
				pdep_candidates[0] = _pdep_u64(number, pdep_mask);
				number = *(input + 1);
				pdep_candidates[1] = _pdep_u64(number, pdep_mask);
				number = *(input + 2);
				pdep_candidates[2] = _pdep_u64(number, pdep_mask);
				number = *(input + 3);
				pdep_candidates[3] = _pdep_u64(number, pdep_mask);

				const uint256_t candidates = _mm256_loadu_si256((uint256_t*)pdep_candidates);

				// select the high and low halves of each byte
				const uint256_t hi_bits = _mm256_and_si256(candidates, nybble_mask);
				const uint256_t lo_bit = _mm256_and_si256(_mm256_srli_epi64(candidates, 4), nybble_mask);

				// replace nybbles with the sum of their remainders
				uint256_t rems_b3 = _mm256_shuffle_epi8(nybble_lookup_b3, hi_bits);
				uint256_t rems_b4 = _mm256_shuffle_epi8(nybble_lookup_b4, hi_bits);
				uint256_t rems_b5 = _mm256_shuffle_epi8(nybble_lookup_b5, hi_bits);
				uint256_t rems_b9 = _mm256_shuffle_epi8(nybble_lookup_b9, hi_bits);

				// add bit 0 to bits 1-4
				rems_b3 = _mm256_add_epi8(rems_b3, lo_bit);
				rems_b4 = _mm256_add_epi8(rems_b4, lo_bit);
				rems_b5 = _mm256_add_epi8(rems_b5, lo_bit);
				rems_b9 = _mm256_add_epi8(rems_b9, lo_bit);
				// h-sum remainders
				rems_b3 = _mm256_sad_epu8(rems_b3, _mm256_setzero_si256());
				rems_b4 = _mm256_sad_epu8(rems_b4, _mm256_setzero_si256());
				rems_b5 = _mm256_sad_epu8(rems_b5, _mm256_setzero_si256());
				rems_b9 = _mm256_sad_epu8(rems_b9, _mm256_setzero_si256());

				const uint256_t rems_b34 = _mm256_packus_epi32(rems_b3, rems_b4);
				const uint256_t rems_b59 = _mm256_packus_epi32(rems_b5, rems_b9);
				const uint256_t rems_b3459 = _mm256_packus_epi32(rems_b34, rems_b59);

				// store results on the stack
				_mm256_storeu_si256((uint256_t*)sums, rems_b3459);
			}

			// run pdep instructions two iterations ahead
			{
				uint64_t number = *(input + 4 + 0);
				pdep_candidates[0] = _pdep_u64(number, pdep_mask);
				number = *(input + 4 + 1);
				pdep_candidates[1] = _pdep_u64(number, pdep_mask);
				number = *(input + 4 + 2);
				pdep_candidates[2] = _pdep_u64(number, pdep_mask);
				number = *(input + 4 + 3);
				pdep_candidates[3] = _pdep_u64(number, pdep_mask);
			}

			// load candidates for the next iteration
			uint256_t candidates = _mm256_loadu_si256((uint256_t*)pdep_candidates);

			for (; input < rounded_end; )
			{
				// load four candidates and run pdep instructions two iterations ahead
				uint64_t number = 0;
				number = *(input + 8 + 0);
				pdep_candidates[0] = _pdep_u64(number, pdep_mask);
				number = *(input + 8 + 1);
				pdep_candidates[1] = _pdep_u64(number, pdep_mask);
				number = *(input + 8 + 2);
				pdep_candidates[2] = _pdep_u64(number, pdep_mask);
				number = *(input + 8 + 3);
				pdep_candidates[3] = _pdep_u64(number, pdep_mask);

				// run vector instructions one iteration ahead

				const uint256_t hi_bits = _mm256_and_si256(candidates, nybble_mask);
				const uint256_t lo_bit = _mm256_and_si256(_mm256_srli_epi64(candidates, 4), nybble_mask);

				// load candidates for the next iteration
				candidates = _mm256_loadu_si256((uint256_t*)pdep_candidates);

				uint256_t rems_b3 = _mm256_shuffle_epi8(nybble_lookup_b3, hi_bits);
				uint256_t rems_b4 = _mm256_shuffle_epi8(nybble_lookup_b4, hi_bits);
				uint256_t rems_b5 = _mm256_shuffle_epi8(nybble_lookup_b5, hi_bits);
				uint256_t rems_b9 = _mm256_shuffle_epi8(nybble_lookup_b9, hi_bits);

				rems_b3 = _mm256_add_epi8(rems_b3, lo_bit);
				rems_b4 = _mm256_add_epi8(rems_b4, lo_bit);
				rems_b5 = _mm256_add_epi8(rems_b5, lo_bit);
				rems_b9 = _mm256_add_epi8(rems_b9, lo_bit);
				rems_b3 = _mm256_sad_epu8(rems_b3, _mm256_setzero_si256());
				rems_b4 = _mm256_sad_epu8(rems_b4, _mm256_setzero_si256());
				rems_b5 = _mm256_sad_epu8(rems_b5, _mm256_setzero_si256());
				rems_b9 = _mm256_sad_epu8(rems_b9, _mm256_setzero_si256());

				const uint256_t rems_b34 = _mm256_packus_epi32(rems_b3, rems_b4);
				const uint256_t rems_b59 = _mm256_packus_epi32(rems_b5, rems_b9);
				const uint256_t rems_b3459 = _mm256_packus_epi32(rems_b34, rems_b59);

				// Only advance the pointer if the number is still a candidate
				const size_t inc_0 = indivisible_b3[sums[0]]
					& indivisible_b4[sums[2]]
					& indivisible_b5[sums[4]]
					& indivisible_b9[sums[6]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_0);

				const size_t inc_1 = indivisible_b3[sums[0 + 1]]
					& indivisible_b4[sums[2 + 1]]
					& indivisible_b5[sums[4 + 1]]
					& indivisible_b9[sums[6 + 1]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_1);

				const size_t inc_2 = indivisible_b3[sums[0 + 8]]
					& indivisible_b4[sums[2 + 8]]
					& indivisible_b5[sums[4 + 8]]
					& indivisible_b9[sums[6 + 8]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_2);

				const size_t inc_3 = indivisible_b3[sums[0 + 1 + 8]]
					& indivisible_b4[sums[2 + 1 + 8]]
					& indivisible_b5[sums[4 + 1 + 8]]
					& indivisible_b9[sums[6 + 1 + 8]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_3);

				// store above results on the stack for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, rems_b3459);
			}

		} // end if (on_fast_path)

		uint64_t number = *input;

		for (; input < candidates_end; )
		{
			const size_t pc_0 = pop_count(number & (bitmask << 0));
			size_t b3_rem = pc_0;
			size_t b4_rem = pc_0;
			size_t b5_rem = pc_0;
			size_t b9_rem = pc_0;
			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b3_rem += pc_1 * pow_mod<3, 1, 11>::rem;
			b4_rem += pc_1 * pow_mod<4, 1, 11>::rem;
			b5_rem += pc_1 * pow_mod<5, 1, 11>::rem;
			b9_rem += pc_1 * pow_mod<9, 1, 11>::rem;
			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b3_rem += pc_2 * pow_mod<3, 2, 11>::rem;
			b4_rem += pc_2 * pow_mod<4, 2, 11>::rem;
			b5_rem += pc_2 * pow_mod<5, 2, 11>::rem;
			b9_rem += pc_2 * pow_mod<9, 2, 11>::rem;
			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b3_rem += pc_3 * pow_mod<3, 3, 11>::rem;
			b4_rem += pc_3 * pow_mod<4, 3, 11>::rem;
			b5_rem += pc_3 * pow_mod<5, 3, 11>::rem;
			b9_rem += pc_3 * pow_mod<9, 3, 11>::rem;
			const size_t pc_4 = pop_count(number & (bitmask << 4));
			b3_rem += pc_4 * pow_mod<3, 4, 11>::rem;
			b4_rem += pc_4 * pow_mod<4, 4, 11>::rem;
			b5_rem += pc_4 * pow_mod<5, 4, 11>::rem;
			b9_rem += pc_4 * pow_mod<9, 4, 11>::rem;

			*output = number; // always write
			number = *++input; // load ahead

			// Only advance the pointer if the number is still a candidate
			size_t inc = indivisible_by_11[b3_rem]
				& indivisible_by_11[b4_rem]
				& indivisible_by_11[b5_rem]
				& indivisible_by_11[b9_rem];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* div_tests_with_10_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// bases 6, 7, and 8 mod 11 (10 remainders)
		constexpr size_t bitmask = bitmask_for<6, 11>::val;
		static_assert(bitmask == bitmask_for<7, 11>::val);
		static_assert(bitmask == bitmask_for<8, 11>::val);
		static_assert(period_of<bitmask>::val == 10);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr size_t upper_bits_mask = size_t(-1) << 32;

			constexpr static uint128_t static_b6_rems_lo = { .m128i_u8{ 1, 6, 3, 7, 9, 10, 5, 8, 4, 2, 1, 6, 3, 7, 9, 10 } };
			constexpr static uint128_t static_b6_rems_hi = { .m128i_u8{ 5, 8, 4, 2, 1, 6, 3, 7, 9, 10, 5, 8, 4, 2, 1, 6 } };
			constexpr static uint128_t static_b7_rems_lo = { .m128i_u8{ 1, 7, 5, 2, 3, 10, 4, 6, 9, 8, 1, 7, 5, 2, 3, 10 } };
			constexpr static uint128_t static_b7_rems_hi = { .m128i_u8{ 4, 6, 9, 8, 1, 7, 5, 2, 3, 10, 4, 6, 9, 8, 1, 7 } };
			constexpr static uint128_t static_b8_rems_lo = { .m128i_u8{ 1, 8, 9, 6, 4, 10, 3, 2, 5, 7, 1, 8, 9, 6, 4, 10 } };
			constexpr static uint128_t static_b8_rems_hi = { .m128i_u8{ 3, 2, 5, 7, 1, 8, 9, 6, 4, 10, 3, 2, 5, 7, 1, 8 } };

			const size_t upper_bits = (*input) & upper_bits_mask;
			const size_t pc_0 = pop_count(upper_bits & (bitmask << 0));
			size_t b6_sum = pc_0;
			size_t b7_sum = pc_0;
			size_t b8_sum = pc_0;
			const size_t pc_1 = pop_count(upper_bits & (bitmask << 1));
			b6_sum += pc_1 * pow_mod<6, 1, 11>::rem;
			b7_sum += pc_1 * pow_mod<7, 1, 11>::rem;
			b8_sum += pc_1 * pow_mod<8, 1, 11>::rem;
			const size_t pc_2 = pop_count(upper_bits & (bitmask << 2));
			b6_sum += pc_2 * pow_mod<6, 2, 11>::rem;
			b7_sum += pc_2 * pow_mod<7, 2, 11>::rem;
			b8_sum += pc_2 * pow_mod<8, 2, 11>::rem;
			const size_t pc_3 = pop_count(upper_bits & (bitmask << 3));
			b6_sum += pc_3 * pow_mod<6, 3, 11>::rem;
			b7_sum += pc_3 * pow_mod<7, 3, 11>::rem;
			b8_sum += pc_3 * pow_mod<8, 3, 11>::rem;
			const size_t pc_4 = pop_count(upper_bits & (bitmask << 4));
			b6_sum += pc_4 * pow_mod<6, 4, 11>::rem;
			b7_sum += pc_4 * pow_mod<7, 4, 11>::rem;
			b8_sum += pc_4 * pow_mod<8, 4, 11>::rem;
			const size_t pc_5 = pop_count(upper_bits & (bitmask << 5));
			b6_sum += pc_5 * pow_mod<6, 5, 11>::rem;
			b7_sum += pc_5 * pow_mod<7, 5, 11>::rem;
			b8_sum += pc_5 * pow_mod<8, 5, 11>::rem;
			const size_t pc_6 = pop_count(upper_bits & (bitmask << 6));
			b6_sum += pc_6 * pow_mod<6, 6, 11>::rem;
			b7_sum += pc_6 * pow_mod<7, 6, 11>::rem;
			b8_sum += pc_6 * pow_mod<8, 6, 11>::rem;
			const size_t pc_7 = pop_count(upper_bits & (bitmask << 7));
			b6_sum += pc_7 * pow_mod<6, 7, 11>::rem;
			b7_sum += pc_7 * pow_mod<7, 7, 11>::rem;
			b8_sum += pc_7 * pow_mod<8, 7, 11>::rem;
			const size_t pc_8 = pop_count(upper_bits & (bitmask << 8));
			b6_sum += pc_8 * pow_mod<6, 8, 11>::rem;
			b7_sum += pc_8 * pow_mod<7, 8, 11>::rem;
			b8_sum += pc_8 * pow_mod<8, 8, 11>::rem;
			const size_t pc_9 = pop_count(upper_bits & (bitmask << 9));
			b6_sum += pc_9 * pow_mod<6, 9, 11>::rem;
			b7_sum += pc_9 * pow_mod<7, 9, 11>::rem;
			b8_sum += pc_9 * pow_mod<8, 9, 11>::rem;

			const uint8_t* const indivisible_b6 = indivisible_by[get_prime_index<11>::idx].data() + b6_sum;
			const uint8_t* const indivisible_b7 = indivisible_by[get_prime_index<11>::idx].data() + b7_sum;
			const uint8_t* const indivisible_b8 = indivisible_by[get_prime_index<11>::idx].data() + b8_sum;

			// calculate the number of candidates, rounded down to the nearest 2
			const size_t n_of_candidates = candidates_end - input;
			const size_t* const rounded_end = input + (n_of_candidates - (n_of_candidates % 2));

			const uint256_t shuffle_mask_lo = _mm256_set_epi64x(0x0909090909090909, 0x0808080808080808, 0x0101010101010101, 0x0000000000000000);
			const uint256_t shuffle_mask_hi = _mm256_set_epi64x(0x0B0B0B0B0B0B0B0B, 0x0A0A0A0A0A0A0A0A, 0x0303030303030303, 0x0202020202020202);
			const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

			const uint128_t xmm_b6_rems_lo = _mm_loadu_si128(&static_b6_rems_lo);
			const uint256_t b6_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b6_rems_lo), xmm_b6_rems_lo, 1);
			const uint128_t xmm_b6_rems_hi = _mm_loadu_si128(&static_b6_rems_hi);
			const uint256_t b6_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b6_rems_hi), xmm_b6_rems_hi, 1);
			const uint128_t xmm_b7_rems_lo = _mm_loadu_si128(&static_b7_rems_lo);
			const uint256_t b7_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b7_rems_lo), xmm_b7_rems_lo, 1);
			const uint128_t xmm_b7_rems_hi = _mm_loadu_si128(&static_b7_rems_hi);
			const uint256_t b7_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b7_rems_hi), xmm_b7_rems_hi, 1);
			const uint128_t xmm_b8_rems_lo = _mm_loadu_si128(&static_b8_rems_lo);
			const uint256_t b8_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b8_rems_lo), xmm_b8_rems_lo, 1);
			const uint128_t xmm_b8_rems_hi = _mm_loadu_si128(&static_b8_rems_hi);
			const uint256_t b8_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b8_rems_hi), xmm_b8_rems_hi, 1);

			uint16_t sums[16]{};

			// run vector instructions one iteration ahead
			{
				// load two candidates
				const uint128_t xmm_candidates = _mm_loadu_si128((uint128_t*)input);
				const uint256_t candidates = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidates), xmm_candidates, 1);

				// convert bits to bytes
				uint256_t candidates_lo = _mm256_shuffle_epi8(candidates, shuffle_mask_lo);
				uint256_t candidates_hi = _mm256_shuffle_epi8(candidates, shuffle_mask_hi);
				candidates_lo = _mm256_andnot_si256(candidates_lo, and_mask);
				candidates_hi = _mm256_andnot_si256(candidates_hi, and_mask);
				candidates_lo = _mm256_cmpeq_epi8(candidates_lo, _mm256_setzero_si256());
				candidates_hi = _mm256_cmpeq_epi8(candidates_hi, _mm256_setzero_si256());

				// select remainders
				const uint256_t b6_sums_lo = _mm256_and_si256(candidates_lo, b6_rems_lo);
				const uint256_t b6_sums_hi = _mm256_and_si256(candidates_hi, b6_rems_hi);
				const uint256_t b7_sums_lo = _mm256_and_si256(candidates_lo, b7_rems_lo);
				const uint256_t b7_sums_hi = _mm256_and_si256(candidates_hi, b7_rems_hi);
				const uint256_t b8_sums_lo = _mm256_and_si256(candidates_lo, b8_rems_lo);
				const uint256_t b8_sums_hi = _mm256_and_si256(candidates_hi, b8_rems_hi);

				// add upper and lower 16 rems
				uint256_t b6_sums = _mm256_add_epi8(b6_sums_lo, b6_sums_hi);
				uint256_t b7_sums = _mm256_add_epi8(b7_sums_lo, b7_sums_hi);
				uint256_t b8_sums = _mm256_add_epi8(b8_sums_lo, b8_sums_hi);

				b6_sums = _mm256_sad_epu8(b6_sums, _mm256_setzero_si256());
				b7_sums = _mm256_sad_epu8(b7_sums, _mm256_setzero_si256());
				b8_sums = _mm256_sad_epu8(b8_sums, _mm256_setzero_si256());

				// pack three vectors down to one
				const uint256_t b67_sums = _mm256_packus_epi32(b6_sums, b7_sums);
				const uint256_t b8x_sums = _mm256_packus_epi32(b8_sums, _mm256_setzero_si256());
				uint256_t b678x_sums = _mm256_packus_epi32(b67_sums, b8x_sums);
				// sums are now stored as [0, 0, 1, 1, 2, 2, x, x][0, 0, 1, 1, 2, 2, x, x]
				b678x_sums = _mm256_add_epi16(b678x_sums, _mm256_srli_si256(b678x_sums, 2));
				// sums are now stored as [0, x, 1, x, 2, x, x, x][0, x, 1, x, 2, x, x, x]

				// store on the stack
				_mm256_storeu_si256((uint256_t*)sums, b678x_sums);
			}

			for (; input < rounded_end; )
			{
				const uint128_t xmm_candidates = _mm_loadu_si128((uint128_t*)(input + 2));
				const uint256_t candidates = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidates), xmm_candidates, 1);

				uint256_t candidates_lo = _mm256_shuffle_epi8(candidates, shuffle_mask_lo);
				uint256_t candidates_hi = _mm256_shuffle_epi8(candidates, shuffle_mask_hi);
				candidates_lo = _mm256_andnot_si256(candidates_lo, and_mask);
				candidates_hi = _mm256_andnot_si256(candidates_hi, and_mask);
				candidates_lo = _mm256_cmpeq_epi8(candidates_lo, _mm256_setzero_si256());
				candidates_hi = _mm256_cmpeq_epi8(candidates_hi, _mm256_setzero_si256());

				const uint256_t b6_sums_lo = _mm256_and_si256(candidates_lo, b6_rems_lo);
				const uint256_t b6_sums_hi = _mm256_and_si256(candidates_hi, b6_rems_hi);
				const uint256_t b7_sums_lo = _mm256_and_si256(candidates_lo, b7_rems_lo);
				const uint256_t b7_sums_hi = _mm256_and_si256(candidates_hi, b7_rems_hi);
				const uint256_t b8_sums_lo = _mm256_and_si256(candidates_lo, b8_rems_lo);
				const uint256_t b8_sums_hi = _mm256_and_si256(candidates_hi, b8_rems_hi);

				uint256_t b6_sums = _mm256_add_epi8(b6_sums_lo, b6_sums_hi);
				uint256_t b7_sums = _mm256_add_epi8(b7_sums_lo, b7_sums_hi);
				uint256_t b8_sums = _mm256_add_epi8(b8_sums_lo, b8_sums_hi);
				b6_sums = _mm256_sad_epu8(b6_sums, _mm256_setzero_si256());
				b7_sums = _mm256_sad_epu8(b7_sums, _mm256_setzero_si256());
				b8_sums = _mm256_sad_epu8(b8_sums, _mm256_setzero_si256());

				const uint256_t b67_sums = _mm256_packus_epi32(b6_sums, b7_sums);
				const uint256_t b8x_sums = _mm256_packus_epi32(b8_sums, _mm256_setzero_si256());
				uint256_t b678x_sums = _mm256_packus_epi32(b67_sums, b8x_sums);
				b678x_sums = _mm256_add_epi16(b678x_sums, _mm256_srli_si256(b678x_sums, 2));

				// Only advance the pointer if the number is still a candidate

				const size_t inc_0 = indivisible_b6[sums[0]]
					& indivisible_b7[sums[2]]
					& indivisible_b8[sums[4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_0);

				const size_t inc_1 = indivisible_b6[sums[0 + 8]]
					& indivisible_b7[sums[2 + 8]]
					& indivisible_b8[sums[4 + 8]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_1);

				// store on the stack for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, b678x_sums);
			} // end for each candidate

			// if there is a remaining (odd) candidate, the loop below runs once and handles it

		} // end if on fast path

		const uint8_t* const indivisible_by_11 = indivisible_by[get_prime_index<11>::idx].data();

		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;
			*output = number; // always write

			const size_t pc_0 = pop_count(number & (bitmask << 0));
			size_t b6_sum = pc_0;
			size_t b7_sum = pc_0;
			size_t b8_sum = pc_0;
			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b6_sum += pc_1 * pow_mod<6, 1, 11>::rem;
			b7_sum += pc_1 * pow_mod<7, 1, 11>::rem;
			b8_sum += pc_1 * pow_mod<8, 1, 11>::rem;
			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b6_sum += pc_2 * pow_mod<6, 2, 11>::rem;
			b7_sum += pc_2 * pow_mod<7, 2, 11>::rem;
			b8_sum += pc_2 * pow_mod<8, 2, 11>::rem;
			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b6_sum += pc_3 * pow_mod<6, 3, 11>::rem;
			b7_sum += pc_3 * pow_mod<7, 3, 11>::rem;
			b8_sum += pc_3 * pow_mod<8, 3, 11>::rem;
			const size_t pc_4 = pop_count(number & (bitmask << 4));
			b6_sum += pc_4 * pow_mod<6, 4, 11>::rem;
			b7_sum += pc_4 * pow_mod<7, 4, 11>::rem;
			b8_sum += pc_4 * pow_mod<8, 4, 11>::rem;
			const size_t pc_5 = pop_count(number & (bitmask << 5));
			b6_sum += pc_5 * pow_mod<6, 5, 11>::rem;
			b7_sum += pc_5 * pow_mod<7, 5, 11>::rem;
			b8_sum += pc_5 * pow_mod<8, 5, 11>::rem;
			const size_t pc_6 = pop_count(number & (bitmask << 6));
			b6_sum += pc_6 * pow_mod<6, 6, 11>::rem;
			b7_sum += pc_6 * pow_mod<7, 6, 11>::rem;
			b8_sum += pc_6 * pow_mod<8, 6, 11>::rem;
			const size_t pc_7 = pop_count(number & (bitmask << 7));
			b6_sum += pc_7 * pow_mod<6, 7, 11>::rem;
			b7_sum += pc_7 * pow_mod<7, 7, 11>::rem;
			b8_sum += pc_7 * pow_mod<8, 7, 11>::rem;
			const size_t pc_8 = pop_count(number & (bitmask << 8));
			b6_sum += pc_8 * pow_mod<6, 8, 11>::rem;
			b7_sum += pc_8 * pow_mod<7, 8, 11>::rem;
			b8_sum += pc_8 * pow_mod<8, 8, 11>::rem;
			const size_t pc_9 = pop_count(number & (bitmask << 9));
			b6_sum += pc_9 * pow_mod<6, 9, 11>::rem;
			b7_sum += pc_9 * pow_mod<7, 9, 11>::rem;
			b8_sum += pc_9 * pow_mod<8, 9, 11>::rem;

			// Only advance the pointer if the number is still a candidate
			size_t inc = indivisible_by_11[b6_sum]
				& indivisible_by_11[b7_sum]
				& indivisible_by_11[b8_sum];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* div_tests_with_12_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// bases 6, 7, and 11 mod 13 (12 remainders)
		constexpr size_t bitmask = bitmask_for<6, 13>::val;
		static_assert(bitmask == bitmask_for<7, 13>::val);
		static_assert(bitmask == bitmask_for<11, 13>::val);
		static_assert(period_of<bitmask>::val == 12);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr size_t upper_bits_mask = size_t(-1) << 32;

			constexpr static uint128_t static_b6_rems_lo = { .m128i_u8{ 1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11, 1, 6, 10, 8 } };
			constexpr static uint128_t static_b6_rems_hi = { .m128i_u8{ 9, 2, 12, 7, 3, 5, 4, 11, 1, 6, 10, 8, 9, 2, 12, 7 } };
			constexpr static uint128_t static_b7_rems_lo = { .m128i_u8{ 1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2, 1, 7, 10, 5 } };
			constexpr static uint128_t static_b7_rems_hi = { .m128i_u8{ 9, 11, 12, 6, 3, 8, 4, 2, 1, 7, 10, 5, 9, 11, 12, 6 } };
			constexpr static uint128_t static_b11_rems_lo = { .m128i_u8{ 1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6, 1, 11, 4, 5 } };
			constexpr static uint128_t static_b11_rems_hi = { .m128i_u8{ 3, 7, 12, 2, 9, 8, 10, 6, 1, 11, 4, 5, 3, 7, 12, 2 } };

			const size_t upper_bits = (*input) & upper_bits_mask;
			const size_t pc_0 = pop_count(upper_bits & (bitmask << 0));
			size_t b6_sum = pc_0;
			size_t b7_sum = pc_0;
			size_t b11_sum = pc_0;
			const size_t pc_1 = pop_count(upper_bits & (bitmask << 1));
			b6_sum += pc_1 * pow_mod<6, 1, 13>::rem;
			b7_sum += pc_1 * pow_mod<7, 1, 13>::rem;
			b11_sum += pc_1 * pow_mod<11, 1, 13>::rem;
			const size_t pc_2 = pop_count(upper_bits & (bitmask << 2));
			b6_sum += pc_2 * pow_mod<6, 2, 13>::rem;
			b7_sum += pc_2 * pow_mod<7, 2, 13>::rem;
			b11_sum += pc_2 * pow_mod<11, 2, 13>::rem;
			const size_t pc_3 = pop_count(upper_bits & (bitmask << 3));
			b6_sum += pc_3 * pow_mod<6, 3, 13>::rem;
			b7_sum += pc_3 * pow_mod<7, 3, 13>::rem;
			b11_sum += pc_3 * pow_mod<11, 3, 13>::rem;
			const size_t pc_4 = pop_count(upper_bits & (bitmask << 4));
			b6_sum += pc_4 * pow_mod<6, 4, 13>::rem;
			b7_sum += pc_4 * pow_mod<7, 4, 13>::rem;
			b11_sum += pc_4 * pow_mod<11, 4, 13>::rem;
			const size_t pc_5 = pop_count(upper_bits & (bitmask << 5));
			b6_sum += pc_5 * pow_mod<6, 5, 13>::rem;
			b7_sum += pc_5 * pow_mod<7, 5, 13>::rem;
			b11_sum += pc_5 * pow_mod<11, 5, 13>::rem;
			const size_t pc_6 = pop_count(upper_bits & (bitmask << 6));
			b6_sum += pc_6 * pow_mod<6, 6, 13>::rem;
			b7_sum += pc_6 * pow_mod<7, 6, 13>::rem;
			b11_sum += pc_6 * pow_mod<11, 6, 13>::rem;
			const size_t pc_7 = pop_count(upper_bits & (bitmask << 7));
			b6_sum += pc_7 * pow_mod<6, 7, 13>::rem;
			b7_sum += pc_7 * pow_mod<7, 7, 13>::rem;
			b11_sum += pc_7 * pow_mod<11, 7, 13>::rem;
			const size_t pc_8 = pop_count(upper_bits & (bitmask << 8));
			b6_sum += pc_8 * pow_mod<6, 8, 13>::rem;
			b7_sum += pc_8 * pow_mod<7, 8, 13>::rem;
			b11_sum += pc_8 * pow_mod<11, 8, 13>::rem;
			const size_t pc_9 = pop_count(upper_bits & (bitmask << 9));
			b6_sum += pc_9 * pow_mod<6, 9, 13>::rem;
			b7_sum += pc_9 * pow_mod<7, 9, 13>::rem;
			b11_sum += pc_9 * pow_mod<11, 9, 13>::rem;
			const size_t pc_10 = pop_count(upper_bits & (bitmask << 10));
			b6_sum += pc_10 * pow_mod<6, 10, 13>::rem;
			b7_sum += pc_10 * pow_mod<7, 10, 13>::rem;
			b11_sum += pc_10 * pow_mod<11, 10, 13>::rem;
			const size_t pc_11 = pop_count(upper_bits & (bitmask << 11));
			b6_sum += pc_11 * pow_mod<6, 11, 13>::rem;
			b7_sum += pc_11 * pow_mod<7, 11, 13>::rem;
			b11_sum += pc_11 * pow_mod<11, 11, 13>::rem;

			const uint8_t* const indivisible_b6 = indivisible_by[get_prime_index<13>::idx].data() + b6_sum;
			const uint8_t* const indivisible_b7 = indivisible_by[get_prime_index<13>::idx].data() + b7_sum;
			const uint8_t* const indivisible_b11 = indivisible_by[get_prime_index<13>::idx].data() + b11_sum;

			const size_t n_of_candidates = candidates_end - input;
			const size_t* const rounded_end = input + (n_of_candidates - (n_of_candidates % 2));

			const uint256_t shuffle_mask_lo = _mm256_set_epi64x(0x0909090909090909, 0x0808080808080808, 0x0101010101010101, 0x0000000000000000);
			const uint256_t shuffle_mask_hi = _mm256_set_epi64x(0x0B0B0B0B0B0B0B0B, 0x0A0A0A0A0A0A0A0A, 0x0303030303030303, 0x0202020202020202);
			const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

			const uint128_t xmm_b6_rems_lo = _mm_loadu_si128(&static_b6_rems_lo);
			const uint256_t b6_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b6_rems_lo), xmm_b6_rems_lo, 1);
			const uint128_t xmm_b6_rems_hi = _mm_loadu_si128(&static_b6_rems_hi);
			const uint256_t b6_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b6_rems_hi), xmm_b6_rems_hi, 1);
			const uint128_t xmm_b7_rems_lo = _mm_loadu_si128(&static_b7_rems_lo);
			const uint256_t b7_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b7_rems_lo), xmm_b7_rems_lo, 1);
			const uint128_t xmm_b7_rems_hi = _mm_loadu_si128(&static_b7_rems_hi);
			const uint256_t b7_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b7_rems_hi), xmm_b7_rems_hi, 1);
			const uint128_t xmm_b11_rems_lo = _mm_loadu_si128(&static_b11_rems_lo);
			const uint256_t b11_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b11_rems_lo), xmm_b11_rems_lo, 1);
			const uint128_t xmm_b11_rems_hi = _mm_loadu_si128(&static_b11_rems_hi);
			const uint256_t b11_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b11_rems_hi), xmm_b11_rems_hi, 1);

			uint16_t sums[16]{};

			// run vector instructions one iteration ahead
			{
				// load two candidates
				const uint128_t xmm_candidates = _mm_loadu_si128((uint128_t*)input);
				const uint256_t candidates = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidates), xmm_candidates, 1);

				// convert bits to bytes
				uint256_t candidates_lo = _mm256_shuffle_epi8(candidates, shuffle_mask_lo);
				uint256_t candidates_hi = _mm256_shuffle_epi8(candidates, shuffle_mask_hi);
				candidates_lo = _mm256_andnot_si256(candidates_lo, and_mask);
				candidates_hi = _mm256_andnot_si256(candidates_hi, and_mask);
				candidates_lo = _mm256_cmpeq_epi8(candidates_lo, _mm256_setzero_si256());
				candidates_hi = _mm256_cmpeq_epi8(candidates_hi, _mm256_setzero_si256());

				// select remainders
				const uint256_t b6_sums_lo = _mm256_and_si256(candidates_lo, b6_rems_lo);
				const uint256_t b6_sums_hi = _mm256_and_si256(candidates_hi, b6_rems_hi);
				const uint256_t b7_sums_lo = _mm256_and_si256(candidates_lo, b7_rems_lo);
				const uint256_t b7_sums_hi = _mm256_and_si256(candidates_hi, b7_rems_hi);
				const uint256_t b11_sums_lo = _mm256_and_si256(candidates_lo, b11_rems_lo);
				const uint256_t b11_sums_hi = _mm256_and_si256(candidates_hi, b11_rems_hi);

				// add upper and lower 16 rems
				uint256_t b6_sums = _mm256_add_epi8(b6_sums_lo, b6_sums_hi);
				uint256_t b7_sums = _mm256_add_epi8(b7_sums_lo, b7_sums_hi);
				uint256_t b11_sums = _mm256_add_epi8(b11_sums_lo, b11_sums_hi);

				b6_sums = _mm256_sad_epu8(b6_sums, _mm256_setzero_si256());
				b7_sums = _mm256_sad_epu8(b7_sums, _mm256_setzero_si256());
				b11_sums = _mm256_sad_epu8(b11_sums, _mm256_setzero_si256());

				// pack three vectors down to one
				const uint256_t b67_sums = _mm256_packus_epi32(b6_sums, b7_sums);
				const uint256_t b11x_sums = _mm256_packus_epi32(b11_sums, _mm256_setzero_si256());
				uint256_t b6_7_11_x_sums = _mm256_packus_epi32(b67_sums, b11x_sums);
				// sums are now stored as [0, 0, 1, 1, 2, 2, x, x][0, 0, 1, 1, 2, 2, x, x]
				b6_7_11_x_sums = _mm256_add_epi16(b6_7_11_x_sums, _mm256_srli_si256(b6_7_11_x_sums, 2));
				// sums are now stored as [0, x, 1, x, 2, x, x, x][0, x, 1, x, 2, x, x, x]

				// store on the stack
				_mm256_storeu_si256((uint256_t*)sums, b6_7_11_x_sums);
			}

			for (; input < rounded_end; )
			{
				const uint128_t xmm_candidates = _mm_loadu_si128((uint128_t*)(input + 2));
				const uint256_t candidates = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidates), xmm_candidates, 1);

				uint256_t candidates_lo = _mm256_shuffle_epi8(candidates, shuffle_mask_lo);
				uint256_t candidates_hi = _mm256_shuffle_epi8(candidates, shuffle_mask_hi);
				candidates_lo = _mm256_andnot_si256(candidates_lo, and_mask);
				candidates_hi = _mm256_andnot_si256(candidates_hi, and_mask);
				candidates_lo = _mm256_cmpeq_epi8(candidates_lo, _mm256_setzero_si256());
				candidates_hi = _mm256_cmpeq_epi8(candidates_hi, _mm256_setzero_si256());

				const uint256_t b6_sums_lo = _mm256_and_si256(candidates_lo, b6_rems_lo);
				const uint256_t b6_sums_hi = _mm256_and_si256(candidates_hi, b6_rems_hi);
				const uint256_t b7_sums_lo = _mm256_and_si256(candidates_lo, b7_rems_lo);
				const uint256_t b7_sums_hi = _mm256_and_si256(candidates_hi, b7_rems_hi);
				const uint256_t b11_sums_lo = _mm256_and_si256(candidates_lo, b11_rems_lo);
				const uint256_t b11_sums_hi = _mm256_and_si256(candidates_hi, b11_rems_hi);

				uint256_t b6_sums = _mm256_add_epi8(b6_sums_lo, b6_sums_hi);
				uint256_t b7_sums = _mm256_add_epi8(b7_sums_lo, b7_sums_hi);
				uint256_t b11_sums = _mm256_add_epi8(b11_sums_lo, b11_sums_hi);
				b6_sums = _mm256_sad_epu8(b6_sums, _mm256_setzero_si256());
				b7_sums = _mm256_sad_epu8(b7_sums, _mm256_setzero_si256());
				b11_sums = _mm256_sad_epu8(b11_sums, _mm256_setzero_si256());

				const uint256_t b67_sums = _mm256_packus_epi32(b6_sums, b7_sums);
				const uint256_t b11x_sums = _mm256_packus_epi32(b11_sums, _mm256_setzero_si256());
				uint256_t b6_7_11_x_sums = _mm256_packus_epi32(b67_sums, b11x_sums);
				b6_7_11_x_sums = _mm256_add_epi16(b6_7_11_x_sums, _mm256_srli_si256(b6_7_11_x_sums, 2));

				// Only advance the pointer if the number is still a candidate

				const size_t inc_0 = indivisible_b6[sums[0]]
					& indivisible_b7[sums[2]]
					& indivisible_b11[sums[4]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_0);

				const size_t inc_1 = indivisible_b6[sums[0 + 8]]
					& indivisible_b7[sums[2 + 8]]
					& indivisible_b11[sums[4 + 8]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_1);

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, b6_7_11_x_sums);
			} // end for each candidate

			// if there is a remaining (odd) candidate, the loop below runs once and handles it

		} // end if on fast path

		const uint8_t* const indivisible_by_13 = indivisible_by[get_prime_index<13>::idx].data();

		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;
			*output = number; // always write

			const size_t pc_0 = pop_count(number & (bitmask << 0));
			size_t b6_sum = pc_0;
			size_t b7_sum = pc_0;
			size_t b11_sum = pc_0;
			const size_t pc_1 = pop_count(number & (bitmask << 1));
			b6_sum += pc_1 * pow_mod<6, 1, 13>::rem;
			b7_sum += pc_1 * pow_mod<7, 1, 13>::rem;
			b11_sum += pc_1 * pow_mod<11, 1, 13>::rem;
			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b6_sum += pc_2 * pow_mod<6, 2, 13>::rem;
			b7_sum += pc_2 * pow_mod<7, 2, 13>::rem;
			b11_sum += pc_2 * pow_mod<11, 2, 13>::rem;
			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b6_sum += pc_3 * pow_mod<6, 3, 13>::rem;
			b7_sum += pc_3 * pow_mod<7, 3, 13>::rem;
			b11_sum += pc_3 * pow_mod<11, 3, 13>::rem;
			const size_t pc_4 = pop_count(number & (bitmask << 4));
			b6_sum += pc_4 * pow_mod<6, 4, 13>::rem;
			b7_sum += pc_4 * pow_mod<7, 4, 13>::rem;
			b11_sum += pc_4 * pow_mod<11, 4, 13>::rem;
			const size_t pc_5 = pop_count(number & (bitmask << 5));
			b6_sum += pc_5 * pow_mod<6, 5, 13>::rem;
			b7_sum += pc_5 * pow_mod<7, 5, 13>::rem;
			b11_sum += pc_5 * pow_mod<11, 5, 13>::rem;
			const size_t pc_6 = pop_count(number & (bitmask << 6));
			b6_sum += pc_6 * pow_mod<6, 6, 13>::rem;
			b7_sum += pc_6 * pow_mod<7, 6, 13>::rem;
			b11_sum += pc_6 * pow_mod<11, 6, 13>::rem;
			const size_t pc_7 = pop_count(number & (bitmask << 7));
			b6_sum += pc_7 * pow_mod<6, 7, 13>::rem;
			b7_sum += pc_7 * pow_mod<7, 7, 13>::rem;
			b11_sum += pc_7 * pow_mod<11, 7, 13>::rem;
			const size_t pc_8 = pop_count(number & (bitmask << 8));
			b6_sum += pc_8 * pow_mod<6, 8, 13>::rem;
			b7_sum += pc_8 * pow_mod<7, 8, 13>::rem;
			b11_sum += pc_8 * pow_mod<11, 8, 13>::rem;
			const size_t pc_9 = pop_count(number & (bitmask << 9));
			b6_sum += pc_9 * pow_mod<6, 9, 13>::rem;
			b7_sum += pc_9 * pow_mod<7, 9, 13>::rem;
			b11_sum += pc_9 * pow_mod<11, 9, 13>::rem;
			const size_t pc_10 = pop_count(number & (bitmask << 10));
			b6_sum += pc_10 * pow_mod<6, 10, 13>::rem;
			b7_sum += pc_10 * pow_mod<7, 10, 13>::rem;
			b11_sum += pc_10 * pow_mod<11, 10, 13>::rem;
			const size_t pc_11 = pop_count(number & (bitmask << 11));
			b6_sum += pc_11 * pow_mod<6, 11, 13>::rem;
			b7_sum += pc_11 * pow_mod<7, 11, 13>::rem;
			b11_sum += pc_11 * pow_mod<11, 11, 13>::rem;

			// Only advance the pointer if the number is still a candidate
			const size_t inc = indivisible_by_13[b6_sum]
				& indivisible_by_13[b7_sum]
				& indivisible_by_13[b11_sum];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* div_tests_with_16_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		constexpr static uint256_t static_rems_01{ .m256i_u8{
			1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6,     // base 3 % 17
			1, 5, 8, 6, 13, 14, 2, 10, 16, 12, 9, 11, 4, 3, 15, 7 } }; // base 5 % 17
		constexpr static uint256_t static_rems_23{ .m256i_u8{
			1, 6, 2, 12, 4, 7, 8, 14, 16, 11, 15, 5, 13, 10, 9, 3,     // base 6 % 17
			1, 7, 15, 3, 4, 11, 9, 12, 16, 10, 2, 14, 13, 6, 8, 5 } }; // base 7 % 17
		constexpr static uint256_t static_rems_45{ .m256i_u8{
			1, 10, 15, 14, 4, 6, 9, 5, 16, 7, 2, 3, 13, 11, 8, 12,     // base 10 % 17
			1, 11, 2, 5, 4, 10, 8, 3, 16, 6, 15, 12, 13, 7, 9, 14 } }; // base 11 % 17
		constexpr static uint256_t static_rems_6x{ .m256i_u8{
			1, 12, 8, 11, 13, 3, 2, 7, 16, 5, 9, 6, 4, 14, 15, 10,     // base 12 % 17
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } }; // padding

		const uint8_t* const indivisible_by_17 = indivisible_by[get_prime_index<17>::idx].data();

		size_t* output = input;
		size_t number = *input;

		size_t upper_bits = number & upper_bits_mask;
		// get popcounts of upper 32 bits, in both lanes, in range 0-2
		uint256_t upper_bytes = util::expand_bits_to_bytes(number >> 32);
		upper_bytes = _mm256_and_si256(upper_bytes, _mm256_set1_epi8(0x01)); // convert 0xFF -> 0x01
		upper_bytes = _mm256_add_epi8(upper_bytes, _mm256_permute2x128_si256(upper_bytes, upper_bytes, 1));

		const uint256_t rems_01 = _mm256_loadu_si256(&static_rems_01);
		const uint256_t rems_23 = _mm256_loadu_si256(&static_rems_23);
		const uint256_t rems_45 = _mm256_loadu_si256(&static_rems_45);
		const uint256_t rems_6x = _mm256_loadu_si256(&static_rems_6x);

		const uint256_t shuffle_mask_lo = _mm256_set_epi64x(0x0101010101010101, 0x0000000000000000, 0x0101010101010101, 0x0000000000000000);
		const uint256_t shuffle_mask_hi = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, 0x0303030303030303, 0x0202020202020202);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		uint64_t sums_0213465x[8]{};

		{
			// this loads more than we need, but VBROADCASTI128 is the cheapest way to load to both lanes
			const uint128_t xmm_candidate = _mm_loadu_si128((uint128_t*)input);
			const uint256_t candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

			// convert bits to bytes
			uint256_t candidate_lo = _mm256_shuffle_epi8(candidate, shuffle_mask_lo);
			uint256_t candidate_hi = _mm256_shuffle_epi8(candidate, shuffle_mask_hi);
			candidate_lo = _mm256_andnot_si256(candidate_lo, and_mask);
			candidate_hi = _mm256_andnot_si256(candidate_hi, and_mask);
			candidate_lo = _mm256_cmpeq_epi8(candidate_lo, _mm256_setzero_si256());
			candidate_hi = _mm256_cmpeq_epi8(candidate_hi, _mm256_setzero_si256());
			candidate_lo = _mm256_and_si256(candidate_lo, _mm256_set1_epi8(0x01));
			candidate_hi = _mm256_and_si256(candidate_hi, _mm256_set1_epi8(0x01));

			// get popcounts of lower 32 bits in range 0-2, in both lanes
			const uint256_t pc_lower = _mm256_add_epi8(candidate_lo, candidate_hi);
			// get popcounts of all 64 bits in range 0-4, in both lanes
			const uint256_t pc = _mm256_add_epi8(pc_lower, upper_bytes);

			// multiply 16+16 8-bit values by 16+16 remainders, storing partially summed results as 8+8 uint16_ts
			const uint256_t sums_01 = _mm256_maddubs_epi16(pc, rems_01);
			const uint256_t sums_23 = _mm256_maddubs_epi16(pc, rems_23);
			const uint256_t sums_45 = _mm256_maddubs_epi16(pc, rems_45);
			const uint256_t sums_6x = _mm256_maddubs_epi16(pc, rems_6x);

			// pack into 4x 8-bit integers
			uint256_t sums_0213 = _mm256_packus_epi16(sums_01, sums_23);
			uint256_t sums_465x = _mm256_packus_epi16(sums_45, sums_6x);

			// h-sum into 4 16-bit integers, 2 in each lane
			sums_0213 = _mm256_sad_epu8(sums_0213, _mm256_setzero_si256());
			sums_465x = _mm256_sad_epu8(sums_465x, _mm256_setzero_si256());

			// store on the stack
			_mm256_storeu_si256((uint256_t*)(sums_0213465x + 0), sums_0213);
			_mm256_storeu_si256((uint256_t*)(sums_0213465x + 4), sums_465x);
		}

		// load ahead
		uint128_t xmm_candidate = _mm_loadu_si128((uint128_t*)(input + 1));
		uint256_t candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

		for (; input < candidates_end; )
		{
			if constexpr (!on_fast_path)
			{
				// Upper 32 bits only change every 4B integers - re-calculate when stale
				const size_t next_n = *(input + 1);
				if ((next_n & upper_bits_mask) != upper_bits)
				{
					upper_bits = next_n & upper_bits_mask;
					upper_bytes = util::expand_bits_to_bytes(next_n >> 32);
					upper_bytes = _mm256_and_si256(upper_bytes, _mm256_set1_epi8(0x01));
					upper_bytes = _mm256_add_epi8(upper_bytes, _mm256_permute2x128_si256(upper_bytes, upper_bytes, 1));
				}
			}

			uint256_t candidate_lo = _mm256_shuffle_epi8(candidate, shuffle_mask_lo);
			uint256_t candidate_hi = _mm256_shuffle_epi8(candidate, shuffle_mask_hi);

			// load two iterations ahead
			xmm_candidate = _mm_loadu_si128((uint128_t*)(input + 2));
			candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

			candidate_lo = _mm256_andnot_si256(candidate_lo, and_mask);
			candidate_hi = _mm256_andnot_si256(candidate_hi, and_mask);
			candidate_lo = _mm256_cmpeq_epi8(candidate_lo, _mm256_setzero_si256());
			candidate_hi = _mm256_cmpeq_epi8(candidate_hi, _mm256_setzero_si256());
			candidate_lo = _mm256_and_si256(candidate_lo, _mm256_set1_epi8(0x01));
			candidate_hi = _mm256_and_si256(candidate_hi, _mm256_set1_epi8(0x01));

			const uint256_t pc_lower = _mm256_add_epi8(candidate_lo, candidate_hi);
			const uint256_t pc = _mm256_add_epi8(pc_lower, upper_bytes);

			const uint256_t sums_01 = _mm256_maddubs_epi16(pc, rems_01);
			const uint256_t sums_23 = _mm256_maddubs_epi16(pc, rems_23);
			const uint256_t sums_45 = _mm256_maddubs_epi16(pc, rems_45);
			const uint256_t sums_6x = _mm256_maddubs_epi16(pc, rems_6x);

			uint256_t sums_0213 = _mm256_packus_epi16(sums_01, sums_23);
			uint256_t sums_465x = _mm256_packus_epi16(sums_45, sums_6x);

			sums_0213 = _mm256_sad_epu8(sums_0213, _mm256_setzero_si256());
			sums_465x = _mm256_sad_epu8(sums_465x, _mm256_setzero_si256());

			*output = *input++; // always write

			// Only advance the pointer if the nth bit is 0 in all lookups
			const size_t inc = indivisible_by_17[sums_0213465x[0]]
				& indivisible_by_17[sums_0213465x[1]]
				& indivisible_by_17[sums_0213465x[2]]
				& indivisible_by_17[sums_0213465x[3]]
				& indivisible_by_17[sums_0213465x[4]]
				& indivisible_by_17[sums_0213465x[5]]
				& indivisible_by_17[sums_0213465x[6]];
			output = (uint64_t*)(((uint8_t*)output) + inc);

			// store the above results for the next iteration
			_mm256_storeu_si256((uint256_t*)(sums_0213465x + 0), sums_0213);
			_mm256_storeu_si256((uint256_t*)(sums_0213465x + 4), sums_465x);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* base13_mod17_div_test(size_t* input,
													const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t bitmask = bitmask_for<13, 17>::val;
		static_assert(period_of<bitmask>::val == 4);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static uint256_t static_nybble_lookup = mbp::detail::build_4rem_shuffle_lookup<13, 17>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup = _mm256_loadu_si256(&static_nybble_lookup);

			const uint8_t* indivisible_lookup_ptr = indivisible_lookup_b13m17.data();
			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			indivisible_lookup_ptr += get_upper_sum_of_rems<17, in_base<13>>(*input);

			alignas(32) uint32_t sums[8] = {};

			// run vector instructions one iteration ahead

			// load 8 candidates into two registers
			uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
			uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));

			// we only want the bottom 32 bits of our 8 candidates
			ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

			// select the high and low halves of each byte
			uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
			uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

			// load 8 candidates into two registers, two iterations ahead
			ymm0 = _mm256_loadu_si256((uint256_t*)(input + 8));
			ymm1 = _mm256_loadu_si256((uint256_t*)(input + 12));

			// replace the bits of each nybble with their remainder
			lo_nybbles = _mm256_shuffle_epi8(nybble_lookup, lo_nybbles);
			hi_nybbles = _mm256_shuffle_epi8(nybble_lookup, hi_nybbles);

			// repack, then horizontally add 8 remainders for each candidate
			uint256_t candidates_0426 = _mm256_unpacklo_epi32(lo_nybbles, hi_nybbles);
			uint256_t candidates_1537 = _mm256_unpackhi_epi32(lo_nybbles, hi_nybbles);
			uint256_t sums_0426 = _mm256_sad_epu8(candidates_0426, _mm256_setzero_si256());
			uint256_t sums_1537 = _mm256_sad_epu8(candidates_1537, _mm256_setzero_si256());

			uint256_t sums_04152637 = _mm256_packus_epi32(sums_0426, sums_1537);

			// store on the stack
			_mm256_storeu_si256((uint256_t*)sums, sums_04152637);

			// break a dependency chain bottleneck by writing to two addresses
			size_t* output_a = output;
			size_t* output_b = output + 1;

			// run vector instructions one iteration ahead; load two iterations ahead
			for (; input < candidates_end_rounded; )
			{
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				ymm0 = _mm256_loadu_si256((uint256_t*)(input + 16));
				ymm1 = _mm256_loadu_si256((uint256_t*)(input + 20));

				lo_nybbles = _mm256_shuffle_epi8(nybble_lookup, lo_nybbles);
				hi_nybbles = _mm256_shuffle_epi8(nybble_lookup, hi_nybbles);

				candidates_0426 = _mm256_unpacklo_epi32(lo_nybbles, hi_nybbles);
				candidates_1537 = _mm256_unpackhi_epi32(lo_nybbles, hi_nybbles);
				sums_0426 = _mm256_sad_epu8(candidates_0426, _mm256_setzero_si256());
				sums_1537 = _mm256_sad_epu8(candidates_1537, _mm256_setzero_si256());

				sums_04152637 = _mm256_packus_epi32(sums_0426, sums_1537);

				const size_t sum_a = sums[0];
				*output_a = *input++; // always copy
				// branchless conditional increment
				output_a = (uint64_t*)(((uint8_t*)output_a) + indivisible_lookup_ptr[sum_a]);

				const size_t sum_b = sums[2];
				*output_b = *input++;
				output_b = (uint64_t*)(((uint8_t*)output_b) + indivisible_lookup_ptr[sum_b]);

				const size_t sum_c = sums[4];
				*output_a = *input++;
				output_a = (uint64_t*)(((uint8_t*)output_a) + indivisible_lookup_ptr[sum_c]);

				const size_t sum_d = sums[6];
				*output_b = *input++;
				output_b = (uint64_t*)(((uint8_t*)output_b) + indivisible_lookup_ptr[sum_d]);

				const size_t sum_e = sums[1];
				*output_a = *input++;
				output_a = (uint64_t*)(((uint8_t*)output_a) + indivisible_lookup_ptr[sum_e]);

				const size_t sum_f = sums[3];
				*output_b = *input++;
				output_b = (uint64_t*)(((uint8_t*)output_b) + indivisible_lookup_ptr[sum_f]);

				const size_t sum_g = sums[5];
				*output_a = *input++;
				output_a = (uint64_t*)(((uint8_t*)output_a) + indivisible_lookup_ptr[sum_g]);

				const size_t sum_h = sums[7];
				*output_b = *input++;
				output_b = (uint64_t*)(((uint8_t*)output_b) + indivisible_lookup_ptr[sum_h]);

				// store on the stack for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, sums_04152637);
			} // end for

			// When the output ptrs end up nonadjacent, move data from a->b or b->a
			// until our data is contiguous again.

			while (output_a > output_b + 1)
			{
				output_a -= 2;
				*output_b = *output_a;
				output_b += 2;
			}

			while (output_a < output_b - 1)
			{
				output_b -= 2;
				*output_a = *output_b;
				output_a += 2;
			}

			// keep the smaller of the two
			output = (output_a < output_b) ? output_a : output_b;

		} // end if on_fast_path

		const uint8_t* const indivisible_lookup_ptr = indivisible_by[get_prime_index<17>::idx].data();

		// handle any elements not handled by the fast loop
		for (; input < candidates_end; ++input)
		{
			// read, then always write
			const size_t candidate = *input;
			*output = candidate;

			size_t sum = 0;
			sum += pop_count(candidate & (bitmask << 0)) * pow_mod<13, 0, 17>::rem;
			sum += pop_count(candidate & (bitmask << 1)) * pow_mod<13, 1, 17>::rem;
			sum += pop_count(candidate & (bitmask << 2)) * pow_mod<13, 2, 17>::rem;
			sum += pop_count(candidate & (bitmask << 3)) * pow_mod<13, 3, 17>::rem;

			// conditionally increment
			output = (uint64_t*)(((uint8_t*)output) + indivisible_lookup_ptr[sum]);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* two_div_tests_with_four_rems(size_t* input,
														   const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t bitmask = bitmask_for<12, 29>::val;
		static_assert(bitmask == bitmask_for<6, 37>::val);
		static_assert(period_of<bitmask>::val == 4);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static uint256_t static_nybble_lookup_12m29 = mbp::detail::build_4rem_shuffle_lookup<12, 29>();
			constexpr static uint256_t static_nybble_lookup_6m37 = mbp::detail::build_4rem_shuffle_lookup<6, 37>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup_12m29 = _mm256_loadu_si256(&static_nybble_lookup_12m29);
			const uint256_t nybble_lookup_6m37 = _mm256_loadu_si256(&static_nybble_lookup_6m37);

			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			const size_t number_upper = (*input) & (size_t(-1) << 32);
			const size_t upc_0 = pop_count(number_upper & (bitmask << 0));
			size_t b12m29_urem = upc_0;
			size_t b6m37_urem = upc_0;
			const size_t upc_1 = pop_count(number_upper & (bitmask << 1));
			b12m29_urem += upc_1 * pow_mod<12, 1, 29>::rem;
			b6m37_urem += upc_1 * pow_mod<6, 1, 37>::rem;
			const size_t upc_2 = pop_count(number_upper & (bitmask << 2));
			b12m29_urem += upc_2 * pow_mod<12, 2, 29>::rem;
			b6m37_urem += upc_2 * pow_mod<6, 2, 37>::rem;
			const size_t upc_3 = pop_count(number_upper & (bitmask << 3));
			b12m29_urem += upc_3 * pow_mod<12, 3, 29>::rem;
			b6m37_urem += upc_3 * pow_mod<6, 3, 37>::rem;
			const uint8_t* const indivisible_12m29 = indivisible_by[get_prime_index<29>::idx].data() + b12m29_urem;
			const uint8_t* const indivisible_6m37 = indivisible_by[get_prime_index<37>::idx].data() + b6m37_urem;

			alignas(32) uint16_t sums[16]{}; // 0426 1537 0426 1537

			// run vector instructions one iteration ahead
			{
				// load 8 candidates into two registers
				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));

				// we only want the bottom 32 bits of our 8 candidates
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				// select the high and low halves of each byte
				const uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				const uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// replace the bits of each nybble with their remainder
				const uint256_t lo_nybbles_12m29 = _mm256_shuffle_epi8(nybble_lookup_12m29, lo_nybbles);
				const uint256_t hi_nybbles_12m29 = _mm256_shuffle_epi8(nybble_lookup_12m29, hi_nybbles);
				const uint256_t lo_nybbles_6m37 = _mm256_shuffle_epi8(nybble_lookup_6m37, lo_nybbles);
				const uint256_t hi_nybbles_6m37 = _mm256_shuffle_epi8(nybble_lookup_6m37, hi_nybbles);
				// repack
				const uint256_t candidates_0426_12m29 = _mm256_unpacklo_epi32(lo_nybbles_12m29, hi_nybbles_12m29);
				const uint256_t candidates_1537_12m29 = _mm256_unpackhi_epi32(lo_nybbles_12m29, hi_nybbles_12m29);
				const uint256_t candidates_0426_6m37 = _mm256_unpacklo_epi32(lo_nybbles_6m37, hi_nybbles_6m37);
				const uint256_t candidates_1537_6m37 = _mm256_unpackhi_epi32(lo_nybbles_6m37, hi_nybbles_6m37);
				// h-sum each 8 remainders -> 1 sum
				const uint256_t sums_0426_12m29 = _mm256_sad_epu8(candidates_0426_12m29, _mm256_setzero_si256());
				const uint256_t sums_1537_12m29 = _mm256_sad_epu8(candidates_1537_12m29, _mm256_setzero_si256());
				const uint256_t sums_0426_6m37 = _mm256_sad_epu8(candidates_0426_6m37, _mm256_setzero_si256());
				const uint256_t sums_1537_6m37 = _mm256_sad_epu8(candidates_1537_6m37, _mm256_setzero_si256());

				uint256_t sums_0404_2626 = _mm256_packus_epi32(sums_0426_12m29, sums_0426_6m37);
				uint256_t sums_1515_3737 = _mm256_packus_epi32(sums_1537_12m29, sums_1537_6m37);
				uint256_t sums_0404_1515_2626_3737 = _mm256_packus_epi32(sums_0404_2626, sums_1515_3737);

				// store results on the stack
				_mm256_storeu_si256((uint256_t*)sums, sums_0404_1515_2626_3737);
			}

			// load two iterations ahead
			uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 8));
			uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 12));

			for (; input < candidates_end_rounded; )
			{
				// run vector instructions one iteration ahead

				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				const uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				const uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// load two iterations ahead
				ymm0 = _mm256_loadu_si256((uint256_t*)(input + 16));
				ymm1 = _mm256_loadu_si256((uint256_t*)(input + 20));

				const uint256_t lo_nybbles_12m29 = _mm256_shuffle_epi8(nybble_lookup_12m29, lo_nybbles);
				const uint256_t hi_nybbles_12m29 = _mm256_shuffle_epi8(nybble_lookup_12m29, hi_nybbles);
				const uint256_t lo_nybbles_6m37 = _mm256_shuffle_epi8(nybble_lookup_6m37, lo_nybbles);
				const uint256_t hi_nybbles_6m37 = _mm256_shuffle_epi8(nybble_lookup_6m37, hi_nybbles);
				const uint256_t candidates_0426_12m29 = _mm256_unpacklo_epi32(lo_nybbles_12m29, hi_nybbles_12m29);
				const uint256_t candidates_1537_12m29 = _mm256_unpackhi_epi32(lo_nybbles_12m29, hi_nybbles_12m29);
				const uint256_t candidates_0426_6m37 = _mm256_unpacklo_epi32(lo_nybbles_6m37, hi_nybbles_6m37);
				const uint256_t candidates_1537_6m37 = _mm256_unpackhi_epi32(lo_nybbles_6m37, hi_nybbles_6m37);
				const uint256_t sums_0426_12m29 = _mm256_sad_epu8(candidates_0426_12m29, _mm256_setzero_si256());
				const uint256_t sums_1537_12m29 = _mm256_sad_epu8(candidates_1537_12m29, _mm256_setzero_si256());
				const uint256_t sums_0426_6m37 = _mm256_sad_epu8(candidates_0426_6m37, _mm256_setzero_si256());
				const uint256_t sums_1537_6m37 = _mm256_sad_epu8(candidates_1537_6m37, _mm256_setzero_si256());

				const uint256_t sums_0404_2626 = _mm256_packus_epi32(sums_0426_12m29, sums_0426_6m37);
				const uint256_t sums_1515_3737 = _mm256_packus_epi32(sums_1537_12m29, sums_1537_6m37);
				const uint256_t sums_0404_1515_2626_3737 = _mm256_packus_epi32(sums_0404_2626, sums_1515_3737);

				*output = *input++; // always copy
				const size_t inc_a =
					indivisible_12m29[sums[0]] &
					indivisible_6m37[sums[0 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_a); // branchless conditional increment

				*output = *input++;
				const size_t inc_b =
					indivisible_12m29[sums[4]] &
					indivisible_6m37[sums[4 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_b);

				*output = *input++;
				const size_t inc_c =
					indivisible_12m29[sums[8]] &
					indivisible_6m37[sums[8 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_c);

				*output = *input++;
				const size_t inc_d =
					indivisible_12m29[sums[12]] &
					indivisible_6m37[sums[12 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_d);

				*output = *input++;
				const size_t inc_e =
					indivisible_12m29[sums[1]] &
					indivisible_6m37[sums[1 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_e);

				*output = *input++;
				const size_t inc_f =
					indivisible_12m29[sums[5]] &
					indivisible_6m37[sums[5 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_f);

				*output = *input++;
				const size_t inc_g =
					indivisible_12m29[sums[9]] &
					indivisible_6m37[sums[9 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_g);

				*output = *input++;
				const size_t inc_h =
					indivisible_12m29[sums[13]] &
					indivisible_6m37[sums[13 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_h);

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, sums_0404_1515_2626_3737);
			} // end for
		} // end if on_fast_path

		const uint8_t* const indivisible_12m29 = indivisible_by[get_prime_index<29>::idx].data();
		const uint8_t* const indivisible_6m37 = indivisible_by[get_prime_index<37>::idx].data();

		// handle any elements not handled by the fast loop
		for (; input < candidates_end; ++input)
		{
			// read, then always write
			const size_t candidate = *input;
			*output = candidate;

			const size_t pc_0 = pop_count(candidate & (bitmask << 0));
			const size_t pc_1 = pop_count(candidate & (bitmask << 1));
			const size_t pc_2 = pop_count(candidate & (bitmask << 2));
			const size_t pc_3 = pop_count(candidate & (bitmask << 3));

			size_t sum_12m29 = pc_0;
			sum_12m29 += pc_1 * pow_mod<12, 1, 29>::rem;
			sum_12m29 += pc_2 * pow_mod<12, 2, 29>::rem;
			sum_12m29 += pc_3 * pow_mod<12, 3, 29>::rem;

			size_t sum_6m37 = pc_0;
			sum_6m37 += pc_1 * pow_mod<6, 1, 37>::rem;
			sum_6m37 += pc_2 * pow_mod<6, 2, 37>::rem;
			sum_6m37 += pc_3 * pow_mod<6, 3, 37>::rem;

			// conditionally increment
			const size_t inc =
				indivisible_12m29[sum_12m29] &
				indivisible_6m37[sum_6m37];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
	tests_are_inlined size_t* div_tests_with_8_rems(size_t* input,
													const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t bitmask = bitmask_for<8, 17>::val;
		static_assert(bitmask == bitmask_for<9, 17>::val);
		static_assert(period_of<bitmask>::val == 8);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static uint256_t static_nybble_lookup_8m17_lo = mbp::detail::build_8rem_shuffle_lookup_lo_nybble<8, 17>();
			constexpr static uint256_t static_nybble_lookup_8m17_hi = mbp::detail::build_8rem_shuffle_lookup_hi_nybble<8, 17>();
			constexpr static uint256_t static_nybble_lookup_9m17_lo = mbp::detail::build_8rem_shuffle_lookup_lo_nybble<9, 17>();
			constexpr static uint256_t static_nybble_lookup_9m17_hi = mbp::detail::build_8rem_shuffle_lookup_hi_nybble<9, 17>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup_8m17_lo = _mm256_loadu_si256(&static_nybble_lookup_8m17_lo);
			const uint256_t nybble_lookup_8m17_hi = _mm256_loadu_si256(&static_nybble_lookup_8m17_hi);
			const uint256_t nybble_lookup_9m17_lo = _mm256_loadu_si256(&static_nybble_lookup_9m17_lo);
			const uint256_t nybble_lookup_9m17_hi = _mm256_loadu_si256(&static_nybble_lookup_9m17_hi);

			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			const size_t number_upper = (*input) & (size_t(-1) << 32);
			const size_t pc_0 = pop_count(number_upper & (bitmask << 0));
			size_t b8m17_rem = pc_0;
			size_t b9m17_rem = pc_0;
			const size_t pc_1 = pop_count(number_upper & (bitmask << 1));
			b8m17_rem += pc_1 * pow_mod<8, 1, 17>::rem;
			b9m17_rem += pc_1 * pow_mod<9, 1, 17>::rem;
			const size_t pc_2 = pop_count(number_upper & (bitmask << 2));
			b8m17_rem += pc_2 * pow_mod<8, 2, 17>::rem;
			b9m17_rem += pc_2 * pow_mod<9, 2, 17>::rem;
			const size_t pc_3 = pop_count(number_upper & (bitmask << 3));
			b8m17_rem += pc_3 * pow_mod<8, 3, 17>::rem;
			b9m17_rem += pc_3 * pow_mod<9, 3, 17>::rem;
			const size_t pc_4 = pop_count(number_upper & (bitmask << 4));
			b8m17_rem += pc_4 * pow_mod<8, 4, 17>::rem;
			b9m17_rem += pc_4 * pow_mod<9, 4, 17>::rem;
			const size_t pc_5 = pop_count(number_upper & (bitmask << 5));
			b8m17_rem += pc_5 * pow_mod<8, 5, 17>::rem;
			b9m17_rem += pc_5 * pow_mod<9, 5, 17>::rem;
			const size_t pc_6 = pop_count(number_upper & (bitmask << 6));
			b8m17_rem += pc_6 * pow_mod<8, 6, 17>::rem;
			b9m17_rem += pc_6 * pow_mod<9, 6, 17>::rem;
			const size_t pc_7 = pop_count(number_upper & (bitmask << 7));
			b8m17_rem += pc_7 * pow_mod<8, 7, 17>::rem;
			b9m17_rem += pc_7 * pow_mod<9, 7, 17>::rem;
			const uint8_t* const indivisible_b8 = indivisible_by[get_prime_index<17>::idx].data() + b8m17_rem;
			const uint8_t* const indivisible_b9 = indivisible_by[get_prime_index<17>::idx].data() + b9m17_rem;

			uint64_t sums[16]{};

			// run vector instructions one iteration ahead
			{
				// load 8 candidates into two registers
				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));

				// we only want the bottom 32 bits of our 8 candidates
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				// select the high and low halves of each byte
				uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// replace the bits of each nybble with their remainder
				uint256_t lo_nybbles_b8m17 = _mm256_shuffle_epi8(nybble_lookup_8m17_lo, lo_nybbles);
				uint256_t hi_nybbles_b8m17 = _mm256_shuffle_epi8(nybble_lookup_8m17_hi, hi_nybbles);
				uint256_t lo_nybbles_b9m17 = _mm256_shuffle_epi8(nybble_lookup_9m17_lo, lo_nybbles);
				uint256_t hi_nybbles_b9m17 = _mm256_shuffle_epi8(nybble_lookup_9m17_hi, hi_nybbles);
				// repack
				uint256_t candidates_0426_b8m17 = _mm256_unpacklo_epi32(lo_nybbles_b8m17, hi_nybbles_b8m17);
				uint256_t candidates_1537_b8m17 = _mm256_unpackhi_epi32(lo_nybbles_b8m17, hi_nybbles_b8m17);
				uint256_t candidates_0426_b9m17 = _mm256_unpacklo_epi32(lo_nybbles_b9m17, hi_nybbles_b9m17);
				uint256_t candidates_1537_b9m17 = _mm256_unpackhi_epi32(lo_nybbles_b9m17, hi_nybbles_b9m17);
				// h-sum each 8 remainders -> 1 sum
				uint256_t sums_0426_b8m17 = _mm256_sad_epu8(candidates_0426_b8m17, _mm256_setzero_si256());
				uint256_t sums_1537_b8m17 = _mm256_sad_epu8(candidates_1537_b8m17, _mm256_setzero_si256());
				uint256_t sums_0426_b9m17 = _mm256_sad_epu8(candidates_0426_b9m17, _mm256_setzero_si256());
				uint256_t sums_1537_b9m17 = _mm256_sad_epu8(candidates_1537_b9m17, _mm256_setzero_si256());

				// store results on the stack
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_0426_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 4), sums_1537_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_0426_b9m17);
				_mm256_storeu_si256((uint256_t*)(sums + 12), sums_1537_b9m17);
			}

			for (; input < candidates_end_rounded; )
			{
				// run vector instructions one iteration ahead

				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 8));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 12));

				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				uint256_t lo_nybbles_b8m17 = _mm256_shuffle_epi8(nybble_lookup_8m17_lo, lo_nybbles);
				uint256_t hi_nybbles_b8m17 = _mm256_shuffle_epi8(nybble_lookup_8m17_hi, hi_nybbles);
				uint256_t lo_nybbles_b9m17 = _mm256_shuffle_epi8(nybble_lookup_9m17_lo, lo_nybbles);
				uint256_t hi_nybbles_b9m17 = _mm256_shuffle_epi8(nybble_lookup_9m17_hi, hi_nybbles);
				uint256_t candidates_0426_b8m17 = _mm256_unpacklo_epi32(lo_nybbles_b8m17, hi_nybbles_b8m17);
				uint256_t candidates_1537_b8m17 = _mm256_unpackhi_epi32(lo_nybbles_b8m17, hi_nybbles_b8m17);
				uint256_t candidates_0426_b9m17 = _mm256_unpacklo_epi32(lo_nybbles_b9m17, hi_nybbles_b9m17);
				uint256_t candidates_1537_b9m17 = _mm256_unpackhi_epi32(lo_nybbles_b9m17, hi_nybbles_b9m17);
				uint256_t sums_0426_b8m17 = _mm256_sad_epu8(candidates_0426_b8m17, _mm256_setzero_si256());
				uint256_t sums_1537_b8m17 = _mm256_sad_epu8(candidates_1537_b8m17, _mm256_setzero_si256());
				uint256_t sums_0426_b9m17 = _mm256_sad_epu8(candidates_0426_b9m17, _mm256_setzero_si256());
				uint256_t sums_1537_b9m17 = _mm256_sad_epu8(candidates_1537_b9m17, _mm256_setzero_si256());

				*output = *input++; // always copy
				const size_t inc_a = indivisible_b8[sums[0]] & indivisible_b9[sums[0 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_a); // branchless conditional increment

				*output = *input++;
				const size_t inc_b = indivisible_b8[sums[4]] & indivisible_b9[sums[4 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_b);

				*output = *input++;
				const size_t inc_c = indivisible_b8[sums[2]] & indivisible_b9[sums[2 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_c);

				*output = *input++;
				const size_t inc_d = indivisible_b8[sums[6]] & indivisible_b9[sums[6 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_d);

				*output = *input++;
				const size_t inc_e = indivisible_b8[sums[1]] & indivisible_b9[sums[1 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_e);

				*output = *input++;
				const size_t inc_f = indivisible_b8[sums[5]] & indivisible_b9[sums[5 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_f);

				*output = *input++;
				const size_t inc_g = indivisible_b8[sums[3]] & indivisible_b9[sums[3 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_g);

				*output = *input++;
				const size_t inc_h = indivisible_b8[sums[7]] & indivisible_b9[sums[7 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_h);

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_0426_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 4), sums_1537_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_0426_b9m17);
				_mm256_storeu_si256((uint256_t*)(sums + 12), sums_1537_b9m17);
			} // end for
		} // end if on_fast_path

		const uint8_t* const indivisible_by_17 = indivisible_by[get_prime_index<17>::idx].data();

		// handle any elements not handled by the fast loop
		for (; input < candidates_end; ++input)
		{
			// read, then always write
			const size_t candidate = *input;
			*output = candidate;

			const size_t pc_0 = pop_count(candidate & (bitmask << 0));
			size_t b8m17_rem = pc_0;
			size_t b9m17_rem = pc_0;
			const size_t pc_1 = pop_count(candidate & (bitmask << 1));
			b8m17_rem += pc_1 * pow_mod<8, 1, 17>::rem;
			b9m17_rem += pc_1 * pow_mod<9, 1, 17>::rem;
			const size_t pc_2 = pop_count(candidate & (bitmask << 2));
			b8m17_rem += pc_2 * pow_mod<8, 2, 17>::rem;
			b9m17_rem += pc_2 * pow_mod<9, 2, 17>::rem;
			const size_t pc_3 = pop_count(candidate & (bitmask << 3));
			b8m17_rem += pc_3 * pow_mod<8, 3, 17>::rem;
			b9m17_rem += pc_3 * pow_mod<9, 3, 17>::rem;
			const size_t pc_4 = pop_count(candidate & (bitmask << 4));
			b8m17_rem += pc_4 * pow_mod<8, 4, 17>::rem;
			b9m17_rem += pc_4 * pow_mod<9, 4, 17>::rem;
			const size_t pc_5 = pop_count(candidate & (bitmask << 5));
			b8m17_rem += pc_5 * pow_mod<8, 5, 17>::rem;
			b9m17_rem += pc_5 * pow_mod<9, 5, 17>::rem;
			const size_t pc_6 = pop_count(candidate & (bitmask << 6));
			b8m17_rem += pc_6 * pow_mod<8, 6, 17>::rem;
			b9m17_rem += pc_6 * pow_mod<9, 6, 17>::rem;
			const size_t pc_7 = pop_count(candidate & (bitmask << 7));
			b8m17_rem += pc_7 * pow_mod<8, 7, 17>::rem;
			b9m17_rem += pc_7 * pow_mod<9, 7, 17>::rem;

			// conditionally increment
			const auto inc = indivisible_by_17[b8m17_rem] & indivisible_by_17[b9m17_rem];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	// currently unused
	template<bool on_fast_path, size_t base, size_t prime>
	tests_are_inlined size_t* custom_4_rem_div_test(size_t* input,
													const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t bitmask = bitmask_for<base, prime>::val;
		static_assert(period_of<bitmask>::val == 4);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static uint256_t static_nybble_lookup = mbp::detail::build_4rem_shuffle_lookup<base, prime>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup = _mm256_loadu_si256(&static_nybble_lookup);

			const uint32_t* pf_lookup_ptr = prime_factor_lookup.data();
			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			pf_lookup_ptr += get_upper_sum_of_rems<prime, in_base<base>>(*input);

			for (; input < candidates_end_rounded; )
			{
				// load 8 candidates into two registers
				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));

				// we only want the bottom 32 bits of our 8 candidates
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				// select the high and low halves of each byte
				uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// replace the bits of each nybble with their remainder
				lo_nybbles = _mm256_shuffle_epi8(nybble_lookup, lo_nybbles);
				hi_nybbles = _mm256_shuffle_epi8(nybble_lookup, hi_nybbles);

				uint256_t candidates_0426 = _mm256_unpacklo_epi32(lo_nybbles, hi_nybbles);
				uint256_t candidates_1537 = _mm256_unpackhi_epi32(lo_nybbles, hi_nybbles);

				uint256_t sums_0426 = _mm256_sad_epu8(candidates_0426, _mm256_setzero_si256());
				uint256_t sums_1537 = _mm256_sad_epu8(candidates_1537, _mm256_setzero_si256());

				// Extract upper lanes
				uint128_t sums_26 = _mm256_extracti128_si256(sums_0426, 1);
				uint128_t sums_37 = _mm256_extracti128_si256(sums_1537, 1);

				*output = *input++; // always copy
				uint32_t mask_a = pf_lookup_ptr[sums_0426.m256i_u64[0]] >> get_prime_index<prime>::idx;
				output += ~mask_a & 1; // branchless conditional increment

				*output = *input++;
				uint32_t mask_b = pf_lookup_ptr[sums_1537.m256i_u64[0]] >> get_prime_index<prime>::idx;
				output += ~mask_b & 1;

				*output = *input++;
				uint32_t mask_c = pf_lookup_ptr[sums_26.m128i_u64[0]] >> get_prime_index<prime>::idx;
				output += ~mask_c & 1;

				*output = *input++;
				uint32_t mask_d = pf_lookup_ptr[sums_37.m128i_u64[0]] >> get_prime_index<prime>::idx;
				output += ~mask_d & 1;

				*output = *input++;
				uint32_t mask_e = pf_lookup_ptr[sums_0426.m256i_u64[1]] >> get_prime_index<prime>::idx;
				output += ~mask_e & 1;

				*output = *input++;
				uint32_t mask_f = pf_lookup_ptr[sums_1537.m256i_u64[1]] >> get_prime_index<prime>::idx;
				output += ~mask_f & 1;

				*output = *input++;
				uint32_t mask_g = pf_lookup_ptr[sums_26.m128i_u64[1]] >> get_prime_index<prime>::idx;
				output += ~mask_g & 1;

				*output = *input++;
				uint32_t mask_h = pf_lookup_ptr[sums_37.m128i_u64[1]] >> get_prime_index<prime>::idx;
				output += ~mask_h & 1;
			} // end for
		} // end if on_fast_path

		// handle any elements not handled by the fast loop
		for (; input < candidates_end; ++input)
		{
			// read, then always write
			const size_t candidate = *input;
			*output = candidate;

			size_t sum = 0;
			sum += pop_count(candidate & (bitmask << 0)) * pow_mod<base, 0, prime>::rem;
			sum += pop_count(candidate & (bitmask << 1)) * pow_mod<base, 1, prime>::rem;
			sum += pop_count(candidate & (bitmask << 2)) * pow_mod<base, 2, prime>::rem;
			sum += pop_count(candidate & (bitmask << 3)) * pow_mod<base, 3, prime>::rem;

			// conditionally increment
			size_t lookup = prime_factor_lookup[sum] >> get_prime_index<prime>::idx;
			output += ~lookup & 0b1;
		}

		return output;
	}

}
