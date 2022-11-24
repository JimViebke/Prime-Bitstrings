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
				if (i & 0x1) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 0, prime>::rem;
				if (i & 0x2) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 1, prime>::rem;
				if (i & 0x4) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 2, prime>::rem;
				if (i & 0x8) lookup.m256i_u8[i] += pow_mod<base, place_value_start + 3, prime>::rem;
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

		const prime_lookup_t* const pf_lookup_ptr = prime_factor_lookup.data();

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
			const prime_lookup_t* const pf_lookup_ptr_a = pf_lookup_ptr + b3m5_urem;
			const prime_lookup_t* const pf_lookup_ptr_b = pf_lookup_ptr + b5m13_urem;
			const prime_lookup_t* const pf_lookup_ptr_c = pf_lookup_ptr + b8m13_urem;
			const prime_lookup_t* const pf_lookup_ptr_d = pf_lookup_ptr + b4m17_urem;

			uint64_t sums_04261537_a[8]{};
			uint64_t sums_04261537_b[8]{};
			uint64_t sums_04261537_c[8]{};
			uint64_t sums_04261537_d[8]{};

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

				// replace the bits of each nybble with their remainder, then sum remainders
				uint256_t rems_lo = _mm256_shuffle_epi8(nybble_lookup_a, nybbles_lo);
				uint256_t rems_hi = _mm256_shuffle_epi8(nybble_lookup_a, nybbles_hi);
				const uint256_t sums_0426_a = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_a = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				rems_lo = _mm256_shuffle_epi8(nybble_lookup_b, nybbles_lo);
				rems_hi = _mm256_shuffle_epi8(nybble_lookup_b, nybbles_hi);
				const uint256_t sums_0426_b = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_b = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				rems_lo = _mm256_shuffle_epi8(nybble_lookup_c, nybbles_lo);
				rems_hi = _mm256_shuffle_epi8(nybble_lookup_c, nybbles_hi);
				const uint256_t sums_0426_c = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_c = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				rems_lo = _mm256_shuffle_epi8(nybble_lookup_d, nybbles_lo);
				rems_hi = _mm256_shuffle_epi8(nybble_lookup_d, nybbles_hi);
				const uint256_t sums_0426_d = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_d = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());

				// store results on the stack
				_mm256_store_si256((uint256_t*)sums_04261537_a, sums_0426_a);
				_mm256_store_si256((uint256_t*)(sums_04261537_a + 4), sums_1537_a);
				_mm256_store_si256((uint256_t*)sums_04261537_b, sums_0426_b);
				_mm256_store_si256((uint256_t*)(sums_04261537_b + 4), sums_1537_b);
				_mm256_store_si256((uint256_t*)sums_04261537_c, sums_0426_c);
				_mm256_store_si256((uint256_t*)(sums_04261537_c + 4), sums_1537_c);
				_mm256_store_si256((uint256_t*)sums_04261537_d, sums_0426_d);
				_mm256_store_si256((uint256_t*)(sums_04261537_d + 4), sums_1537_d);
			}

			for (; input < candidates_end_rounded; )
			{
				// run vector instructions one iteration ahead

				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 8));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 12));

				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				const uint256_t nybbles_lo = _mm256_and_si256(ymm0, nybble_mask);
				const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				uint256_t rems_lo = _mm256_shuffle_epi8(nybble_lookup_a, nybbles_lo);
				uint256_t rems_hi = _mm256_shuffle_epi8(nybble_lookup_a, nybbles_hi);
				const uint256_t sums_0426_a = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_a = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				rems_lo = _mm256_shuffle_epi8(nybble_lookup_b, nybbles_lo);
				rems_hi = _mm256_shuffle_epi8(nybble_lookup_b, nybbles_hi);
				const uint256_t sums_0426_b = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_b = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				rems_lo = _mm256_shuffle_epi8(nybble_lookup_c, nybbles_lo);
				rems_hi = _mm256_shuffle_epi8(nybble_lookup_c, nybbles_hi);
				const uint256_t sums_0426_c = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_c = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				rems_lo = _mm256_shuffle_epi8(nybble_lookup_d, nybbles_lo);
				rems_hi = _mm256_shuffle_epi8(nybble_lookup_d, nybbles_hi);
				const uint256_t sums_0426_d = _mm256_sad_epu8(_mm256_unpacklo_epi32(rems_lo, rems_hi), _mm256_setzero_si256());
				const uint256_t sums_1537_d = _mm256_sad_epu8(_mm256_unpackhi_epi32(rems_lo, rems_hi), _mm256_setzero_si256());

				// Only advance the pointer if the number is still a candidate, that is,
				// the relevant bit from each lookup is 0

				size_t merged_masks_0 = 0;
				merged_masks_0 |= (pf_lookup_ptr_a[sums_04261537_a[0]] >> get_prime_index<5>::idx);
				merged_masks_0 |= (pf_lookup_ptr_b[sums_04261537_b[0]] >> get_prime_index<13>::idx);
				merged_masks_0 |= (pf_lookup_ptr_c[sums_04261537_c[0]] >> get_prime_index<13>::idx);
				merged_masks_0 |= (pf_lookup_ptr_d[sums_04261537_d[0]] >> get_prime_index<17>::idx);
				*output = *input++; // always copy
				output += ~merged_masks_0 & 0b1; // branchless conditional increment

				size_t merged_masks_1 = 0;
				merged_masks_1 |= (pf_lookup_ptr_a[sums_04261537_a[4]] >> get_prime_index<5>::idx);
				merged_masks_1 |= (pf_lookup_ptr_b[sums_04261537_b[4]] >> get_prime_index<13>::idx);
				merged_masks_1 |= (pf_lookup_ptr_c[sums_04261537_c[4]] >> get_prime_index<13>::idx);
				merged_masks_1 |= (pf_lookup_ptr_d[sums_04261537_d[4]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_1 & 0b1;

				size_t merged_masks_2 = 0;
				merged_masks_2 |= (pf_lookup_ptr_a[sums_04261537_a[2]] >> get_prime_index<5>::idx);
				merged_masks_2 |= (pf_lookup_ptr_b[sums_04261537_b[2]] >> get_prime_index<13>::idx);
				merged_masks_2 |= (pf_lookup_ptr_c[sums_04261537_c[2]] >> get_prime_index<13>::idx);
				merged_masks_2 |= (pf_lookup_ptr_d[sums_04261537_d[2]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_2 & 0b1;

				size_t merged_masks_3 = 0;
				merged_masks_3 |= (pf_lookup_ptr_a[sums_04261537_a[6]] >> get_prime_index<5>::idx);
				merged_masks_3 |= (pf_lookup_ptr_b[sums_04261537_b[6]] >> get_prime_index<13>::idx);
				merged_masks_3 |= (pf_lookup_ptr_c[sums_04261537_c[6]] >> get_prime_index<13>::idx);
				merged_masks_3 |= (pf_lookup_ptr_d[sums_04261537_d[6]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_3 & 0b1;

				size_t merged_masks_4 = 0;
				merged_masks_4 |= (pf_lookup_ptr_a[sums_04261537_a[1]] >> get_prime_index<5>::idx);
				merged_masks_4 |= (pf_lookup_ptr_b[sums_04261537_b[1]] >> get_prime_index<13>::idx);
				merged_masks_4 |= (pf_lookup_ptr_c[sums_04261537_c[1]] >> get_prime_index<13>::idx);
				merged_masks_4 |= (pf_lookup_ptr_d[sums_04261537_d[1]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_4 & 0b1;

				size_t merged_masks_5 = 0;
				merged_masks_5 |= (pf_lookup_ptr_a[sums_04261537_a[5]] >> get_prime_index<5>::idx);
				merged_masks_5 |= (pf_lookup_ptr_b[sums_04261537_b[5]] >> get_prime_index<13>::idx);
				merged_masks_5 |= (pf_lookup_ptr_c[sums_04261537_c[5]] >> get_prime_index<13>::idx);
				merged_masks_5 |= (pf_lookup_ptr_d[sums_04261537_d[5]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_5 & 0b1;

				size_t merged_masks_6 = 0;
				merged_masks_6 |= (pf_lookup_ptr_a[sums_04261537_a[3]] >> get_prime_index<5>::idx);
				merged_masks_6 |= (pf_lookup_ptr_b[sums_04261537_b[3]] >> get_prime_index<13>::idx);
				merged_masks_6 |= (pf_lookup_ptr_c[sums_04261537_c[3]] >> get_prime_index<13>::idx);
				merged_masks_6 |= (pf_lookup_ptr_d[sums_04261537_d[3]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_6 & 0b1;

				size_t merged_masks_7 = 0;
				merged_masks_7 |= (pf_lookup_ptr_a[sums_04261537_a[7]] >> get_prime_index<5>::idx);
				merged_masks_7 |= (pf_lookup_ptr_b[sums_04261537_b[7]] >> get_prime_index<13>::idx);
				merged_masks_7 |= (pf_lookup_ptr_c[sums_04261537_c[7]] >> get_prime_index<13>::idx);
				merged_masks_7 |= (pf_lookup_ptr_d[sums_04261537_d[7]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_7 & 0b1;

				// store on the stack for the next iteration
				_mm256_store_si256((uint256_t*)sums_04261537_a, sums_0426_a);
				_mm256_store_si256((uint256_t*)(sums_04261537_a + 4), sums_1537_a);
				_mm256_store_si256((uint256_t*)sums_04261537_b, sums_0426_b);
				_mm256_store_si256((uint256_t*)(sums_04261537_b + 4), sums_1537_b);
				_mm256_store_si256((uint256_t*)sums_04261537_c, sums_0426_c);
				_mm256_store_si256((uint256_t*)(sums_04261537_c + 4), sums_1537_c);
				_mm256_store_si256((uint256_t*)sums_04261537_d, sums_0426_d);
				_mm256_store_si256((uint256_t*)(sums_04261537_d + 4), sums_1537_d);
			}
		}

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

			size_t merged_masks = 0;
			merged_masks |= (pf_lookup_ptr[b3m5_rem] >> get_prime_index<5>::idx);
			merged_masks |= (pf_lookup_ptr[b5m13_rem] >> get_prime_index<13>::idx);
			merged_masks |= (pf_lookup_ptr[b8m13_rem] >> get_prime_index<13>::idx);
			merged_masks |= (pf_lookup_ptr[b4m17_rem] >> get_prime_index<17>::idx);
			output += ~merged_masks & 0b1;
		}

		return output;
	}

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

		const prime_lookup_t* const prime_factor_lookup_ptr = prime_factor_lookup.data();

		size_t* output = input;

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
			size_t merged_masks = 0;
			merged_masks |= (prime_factor_lookup_ptr[b4_m7_rem] >> get_prime_index<7>::idx);
			merged_masks |= (prime_factor_lookup_ptr[b3_m13_rem] >> get_prime_index<13>::idx);
			merged_masks |= (prime_factor_lookup_ptr[b9_m13_rem] >> get_prime_index<13>::idx);

			output += ~merged_masks & 0b1;
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

		constexpr static uint256_t static_rems_01_lo = uint256_t{ .m256i_u8{
			1, 3, 2, 6, 4, 5, 1, 3, 2, 6, 4, 5, 1, 3, 2, 6,
			1, 5, 4, 6, 2, 3, 1, 5, 4, 6, 2, 3, 1, 5, 4, 6 } };
		constexpr static uint256_t static_rems_01_hi = uint256_t{ .m256i_u8{
			4, 5, 1, 3, 2, 6, 4, 5, 1, 3, 2, 6, 4, 5, 1, 3,
			2, 3, 1, 5, 4, 6, 2, 3, 1, 5, 4, 6, 2, 3, 1, 5 } };
		constexpr static uint256_t static_rems_23_lo = uint256_t{ .m256i_u8{
			1, 4, 3, 12, 9, 10, 1, 4, 3, 12, 9, 10, 1, 4, 3, 12,
			1, 10, 9, 12, 3, 4, 1, 10, 9, 12, 3, 4, 1, 10, 9, 12 } };
		constexpr static uint256_t static_rems_23_hi = uint256_t{ .m256i_u8{
			9, 10, 1, 4, 3, 12, 9, 10, 1, 4, 3, 12, 9, 10, 1, 4,
			3, 4, 1, 10, 9, 12, 3, 4, 1, 10, 9, 12, 3, 4, 1, 10 } };

		constexpr static uint256_t static_shuffle_mask_lo = uint256_t{ .m256i_u64{
			0x0000000000000000, 0x0101010101010101, 0x0000000000000000, 0x0101010101010101 } };
		constexpr static uint256_t static_shuffle_mask_hi = uint256_t{ .m256i_u64{
			0x0202020202020202, 0x0303030303030303, 0x0202020202020202, 0x0303030303030303 } };

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

		const prime_lookup_t* pf_lookup_ptr_b3m7 = prime_factor_lookup.data() + upper_sum_0;
		const prime_lookup_t* pf_lookup_ptr_b5m7 = prime_factor_lookup.data() + upper_sum_1;
		const prime_lookup_t* pf_lookup_ptr_b4m13 = prime_factor_lookup.data() + upper_sum_2;
		const prime_lookup_t* pf_lookup_ptr_b10m13 = prime_factor_lookup.data() + upper_sum_3;

		const uint256_t shuffle_mask_lo = _mm256_loadu_si256(&static_shuffle_mask_lo);
		const uint256_t shuffle_mask_hi = _mm256_loadu_si256(&static_shuffle_mask_hi);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		const uint256_t rems_01_lo = _mm256_loadu_si256(&static_rems_01_lo);
		const uint256_t rems_01_hi = _mm256_loadu_si256(&static_rems_01_hi);
		const uint256_t rems_23_lo = _mm256_loadu_si256(&static_rems_23_lo);
		const uint256_t rems_23_hi = _mm256_loadu_si256(&static_rems_23_hi);

		uint32_t sums_0x2x1x3x[8]{};

		// run vector instructions one iteration ahead
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

			// mask to select remainders
			const uint256_t sums_01_lo = _mm256_and_si256(candidate_lo, rems_01_lo);
			const uint256_t sums_01_hi = _mm256_and_si256(candidate_hi, rems_01_hi);
			const uint256_t sums_23_lo = _mm256_and_si256(candidate_lo, rems_23_lo);
			const uint256_t sums_23_hi = _mm256_and_si256(candidate_hi, rems_23_hi);

			// v-sum upper and lower 16 elements
			uint256_t sums_01 = _mm256_add_epi8(sums_01_lo, sums_01_hi);
			uint256_t sums_23 = _mm256_add_epi8(sums_23_lo, sums_23_hi);
			// h-sum 16 partial sums per lane -> 2 partial sums per lane
			sums_01 = _mm256_sad_epu8(sums_01, _mm256_setzero_si256());
			sums_23 = _mm256_sad_epu8(sums_23, _mm256_setzero_si256());
			// pack 16 bits per u64 -> 16 bits per u32
			uint256_t sums_0213 = _mm256_packus_epi32(sums_01, sums_23);
			// shift and add to get [0, x, 2, x, 1, x, 3, x]
			sums_0213 = _mm256_add_epi32(sums_0213, _mm256_srli_si256(sums_0213, 4));

			// store on the stack
			_mm256_storeu_si256((uint256_t*)sums_0x2x1x3x, sums_0213);
		}

		for (; input < candidates_end; )
		{
			if constexpr (!on_fast_path)
			{
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;

					// recalculate
					pf_lookup_ptr_b3m7 = prime_factor_lookup.data() + get_upper_sum_of_rems<7, in_base<3>>(number);
					pf_lookup_ptr_b5m7 = prime_factor_lookup.data() + get_upper_sum_of_rems<7, in_base<5>>(number);
					pf_lookup_ptr_b4m13 = prime_factor_lookup.data() + get_upper_sum_of_rems<13, in_base<4>>(number);
					pf_lookup_ptr_b10m13 = prime_factor_lookup.data() + get_upper_sum_of_rems<13, in_base<10>>(number);
				}
			}

			// this loads more than we need, but VBROADCASTI128 is the cheapest way to load to both lanes
			const uint128_t xmm_candidate = _mm_loadu_si128((uint128_t*)(input + 1));
			const uint256_t candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

			uint256_t candidate_lo = _mm256_shuffle_epi8(candidate, shuffle_mask_lo);
			uint256_t candidate_hi = _mm256_shuffle_epi8(candidate, shuffle_mask_hi);
			candidate_lo = _mm256_andnot_si256(candidate_lo, and_mask);
			candidate_hi = _mm256_andnot_si256(candidate_hi, and_mask);
			candidate_lo = _mm256_cmpeq_epi8(candidate_lo, _mm256_setzero_si256());
			candidate_hi = _mm256_cmpeq_epi8(candidate_hi, _mm256_setzero_si256());

			const uint256_t sums_01_lo = _mm256_and_si256(candidate_lo, rems_01_lo);
			const uint256_t sums_01_hi = _mm256_and_si256(candidate_hi, rems_01_hi);
			const uint256_t sums_23_lo = _mm256_and_si256(candidate_lo, rems_23_lo);
			const uint256_t sums_23_hi = _mm256_and_si256(candidate_hi, rems_23_hi);

			uint256_t sums_01 = _mm256_add_epi8(sums_01_lo, sums_01_hi);
			uint256_t sums_23 = _mm256_add_epi8(sums_23_lo, sums_23_hi);
			sums_01 = _mm256_sad_epu8(sums_01, _mm256_setzero_si256());
			sums_23 = _mm256_sad_epu8(sums_23, _mm256_setzero_si256());
			uint256_t sums_0213 = _mm256_packus_epi32(sums_01, sums_23);
			sums_0213 = _mm256_add_epi32(sums_0213, _mm256_srli_si256(sums_0213, 4));

			*output = number; // always write
			number = *++input; // load one iteration ahead

			// Only advance the pointer if the number is still a candidate
			size_t merged_masks = 0;
			merged_masks |= (pf_lookup_ptr_b3m7[sums_0x2x1x3x[0]] >> get_prime_index<7>::idx);
			merged_masks |= (pf_lookup_ptr_b5m7[sums_0x2x1x3x[4]] >> get_prime_index<7>::idx);
			merged_masks |= (pf_lookup_ptr_b4m13[sums_0x2x1x3x[2]] >> get_prime_index<13>::idx);
			merged_masks |= (pf_lookup_ptr_b10m13[sums_0x2x1x3x[6]] >> get_prime_index<13>::idx);

			output += ~merged_masks & 0b1;

			// store on the stack for the next iteration
			_mm256_storeu_si256((uint256_t*)sums_0x2x1x3x, sums_0213);
		}

		return output;
	}

	tests_are_inlined size_t* div_tests_with_five_rems(size_t* input,
													   const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// Intellisense may generate a number of false positives here
		constexpr size_t bitmask = bitmask_for<3, 11>::val;  // base 3 % 11   (5 remainders)
		static_assert(bitmask == bitmask_for<4, 11>::val);  // base 4 % 11   (5 remainders)
		static_assert(bitmask == bitmask_for<5, 11>::val);  // base 5 % 11   (5 remainders)
		static_assert(bitmask == bitmask_for<9, 11>::val);  // base 9 % 11   (5 remainders)
		static_assert(period_of<bitmask>::val == 5);

		// constexpr uint128_t static_rems0 = { always all ones };
		constexpr static uint128_t static_rems1 = uint128_t{ .m128i_u32{
			pow_mod<3, 1, 11>::rem,
			pow_mod<4, 1, 11>::rem,
			pow_mod<5, 1, 11>::rem,
			pow_mod<9, 1, 11>::rem } };
		constexpr static uint128_t static_rems2 = uint128_t{ .m128i_u32{
			pow_mod<3, 2, 11>::rem,
			pow_mod<4, 2, 11>::rem,
			pow_mod<5, 2, 11>::rem,
			pow_mod<9, 2, 11>::rem } };
		constexpr static uint128_t static_rems3 = uint128_t{ .m128i_u32{
			pow_mod<3, 3, 11>::rem,
			pow_mod<4, 3, 11>::rem,
			pow_mod<5, 3, 11>::rem,
			pow_mod<9, 3, 11>::rem } };
		constexpr static uint128_t static_rems4 = uint128_t{ .m128i_u32{
			pow_mod<3, 4, 11>::rem,
			pow_mod<4, 4, 11>::rem,
			pow_mod<5, 4, 11>::rem,
			pow_mod<9, 4, 11>::rem } };

		const auto xmm_rems1 = _mm_loadu_si128(&static_rems1);
		const auto xmm_rems2 = _mm_loadu_si128(&static_rems2);
		const auto xmm_rems3 = _mm_loadu_si128(&static_rems3);
		const auto xmm_rems4 = _mm_loadu_si128(&static_rems4);

		const prime_lookup_t* const prime_factor_lookup_ptr = prime_factor_lookup.data();

		size_t* output = input;

		uint32_t sums[4]{};

		// run vector instructions one iteration ahead
		{
			const size_t number = *input;
			auto xmm0 = _mm_set1_epi32((int)pop_count(number & (bitmask << 0)));
			auto xmm1 = _mm_set1_epi32((int)pop_count(number & (bitmask << 1)));
			auto xmm2 = _mm_set1_epi32((int)pop_count(number & (bitmask << 2)));
			auto xmm3 = _mm_set1_epi32((int)pop_count(number & (bitmask << 3)));
			auto xmm4 = _mm_set1_epi32((int)pop_count(number & (bitmask << 4)));

			// multiply each popcount with Nth remainder of each test
			// (rems0 is always all ones)
			xmm1 = _mm_mullo_epi32(xmm1, xmm_rems1);
			xmm2 = _mm_mullo_epi32(xmm2, xmm_rems2);
			xmm3 = _mm_mullo_epi32(xmm3, xmm_rems3);
			xmm4 = _mm_mullo_epi32(xmm4, xmm_rems4);

			// add pairs
			xmm0 = _mm_add_epi32(xmm0, xmm1);
			xmm2 = _mm_add_epi32(xmm2, xmm3);
			// final adds
			xmm0 = _mm_add_epi32(xmm0, xmm2);
			xmm0 = _mm_add_epi32(xmm0, xmm4);

			// store the above results for the next iteraton
			_mm_storeu_si128((uint128_t*)sums, xmm0);
		}

		for (; input < candidates_end; )
		{
			// load one iteration ahead
			const size_t next_n = *(input + 1);
			auto xmm0 = _mm_set1_epi32((int)pop_count(next_n & (bitmask << 0)));
			auto xmm1 = _mm_set1_epi32((int)pop_count(next_n & (bitmask << 1)));
			auto xmm2 = _mm_set1_epi32((int)pop_count(next_n & (bitmask << 2)));
			auto xmm3 = _mm_set1_epi32((int)pop_count(next_n & (bitmask << 3)));
			auto xmm4 = _mm_set1_epi32((int)pop_count(next_n & (bitmask << 4)));

			xmm1 = _mm_mullo_epi32(xmm1, xmm_rems1);
			xmm2 = _mm_mullo_epi32(xmm2, xmm_rems2);
			xmm3 = _mm_mullo_epi32(xmm3, xmm_rems3);
			xmm4 = _mm_mullo_epi32(xmm4, xmm_rems4);

			xmm0 = _mm_add_epi32(xmm0, xmm1);
			xmm2 = _mm_add_epi32(xmm2, xmm3);
			xmm0 = _mm_add_epi32(xmm0, xmm2);
			xmm0 = _mm_add_epi32(xmm0, xmm4);

			// Only advance the pointer if the number is still a candidate
			size_t merged_lookups = prime_factor_lookup_ptr[sums[0]];
			merged_lookups |= prime_factor_lookup_ptr[sums[1]];
			merged_lookups |= prime_factor_lookup_ptr[sums[2]];
			merged_lookups |= prime_factor_lookup_ptr[sums[3]];

			merged_lookups >>= get_prime_index<11>::idx;

			*output = *input++;
			output += ~merged_lookups & 0b1;

			// store the above results for the next iteraton
			_mm_storeu_si128((uint128_t*)sums, xmm0);
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

			const prime_lookup_t* pf_lookup_ptr_b6 = prime_factor_lookup.data() + b6_sum;
			const prime_lookup_t* pf_lookup_ptr_b7 = prime_factor_lookup.data() + b7_sum;
			const prime_lookup_t* pf_lookup_ptr_b8 = prime_factor_lookup.data() + b8_sum;

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

				size_t merged_lookups_0 = 0;
				merged_lookups_0 |= pf_lookup_ptr_b6[sums[0]];
				merged_lookups_0 |= pf_lookup_ptr_b7[sums[2]];
				merged_lookups_0 |= pf_lookup_ptr_b8[sums[4]];
				merged_lookups_0 >>= get_prime_index<11>::idx;
				*output = *input++;
				output += ~merged_lookups_0 & 0b1;

				size_t merged_lookups_1 = 0;
				merged_lookups_1 |= pf_lookup_ptr_b6[sums[0 + 8]];
				merged_lookups_1 |= pf_lookup_ptr_b7[sums[2 + 8]];
				merged_lookups_1 |= pf_lookup_ptr_b8[sums[4 + 8]];
				merged_lookups_1 >>= get_prime_index<11>::idx;
				*output = *input++;
				output += ~merged_lookups_1 & 0b1;

				// store on the stack for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, b678x_sums);
			} // end for each candidate

			// if there is a remaining (odd) candidate, the loop below runs once and handles it

		} // end if on fast path

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

			size_t merged_lookups = 0;
			merged_lookups |= prime_factor_lookup[b6_sum];
			merged_lookups |= prime_factor_lookup[b7_sum];
			merged_lookups |= prime_factor_lookup[b8_sum];
			merged_lookups >>= get_prime_index<11>::idx;

			// Only advance the pointer if the number is still a candidate
			output += ~merged_lookups & 0b1;
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

			const prime_lookup_t* pf_lookup_ptr_b6 = prime_factor_lookup.data() + b6_sum;
			const prime_lookup_t* pf_lookup_ptr_b7 = prime_factor_lookup.data() + b7_sum;
			const prime_lookup_t* pf_lookup_ptr_b11 = prime_factor_lookup.data() + b11_sum;

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

				size_t merged_lookups_0 = 0;
				merged_lookups_0 |= pf_lookup_ptr_b6[sums[0]];
				merged_lookups_0 |= pf_lookup_ptr_b7[sums[2]];
				merged_lookups_0 |= pf_lookup_ptr_b11[sums[4]];
				merged_lookups_0 >>= get_prime_index<13>::idx;
				*output = *input++;
				output += ~merged_lookups_0 & 0b1;

				size_t merged_lookups_1 = 0;
				merged_lookups_1 |= pf_lookup_ptr_b6[sums[0 + 8]];
				merged_lookups_1 |= pf_lookup_ptr_b7[sums[2 + 8]];
				merged_lookups_1 |= pf_lookup_ptr_b11[sums[4 + 8]];
				merged_lookups_1 >>= get_prime_index<13>::idx;
				*output = *input++;
				output += ~merged_lookups_1 & 0b1;

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)sums, b6_7_11_x_sums);
			} // end for each candidate

			// if there is a remaining (odd) candidate, the loop below runs once and handles it

		} // end if on fast path

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

			size_t merged_lookups = 0;
			merged_lookups |= prime_factor_lookup[b6_sum];
			merged_lookups |= prime_factor_lookup[b7_sum];
			merged_lookups |= prime_factor_lookup[b11_sum];
			merged_lookups >>= get_prime_index<13>::idx;

			// Only advance the pointer if the number is still a candidate
			output += ~merged_lookups & 0b1;
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

		const prime_lookup_t* const pf_lookup_ptr = prime_factor_lookup.data();

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

			// run vector instructions one iteration ahead

			const uint128_t xmm_candidate = _mm_loadu_si128((uint128_t*)(input + 1));
			const uint256_t candidate = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_candidate), xmm_candidate, 1);

			uint256_t candidate_lo = _mm256_shuffle_epi8(candidate, shuffle_mask_lo);
			uint256_t candidate_hi = _mm256_shuffle_epi8(candidate, shuffle_mask_hi);
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
			prime_lookup_t merged_lookups = 0;
			merged_lookups |= pf_lookup_ptr[sums_0213465x[0]];
			merged_lookups |= pf_lookup_ptr[sums_0213465x[1]];
			merged_lookups |= pf_lookup_ptr[sums_0213465x[2]];
			merged_lookups |= pf_lookup_ptr[sums_0213465x[3]];
			merged_lookups |= pf_lookup_ptr[sums_0213465x[4]];
			merged_lookups |= pf_lookup_ptr[sums_0213465x[5]];
			merged_lookups |= pf_lookup_ptr[sums_0213465x[6]];

			// Invert the bit we care about, move it to the 0th position, and mask to select it
			output += (~merged_lookups >> get_prime_index<17>::idx) & 0b1;

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

			const uint32_t* pf_lookup_ptr = prime_factor_lookup.data();
			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			pf_lookup_ptr += get_upper_sum_of_rems<17, in_base<13>>(*input);

			uint64_t sums_04261537[8] = {};

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

			// store on the stack
			_mm256_storeu_si256((uint256_t*)sums_04261537, sums_0426);
			_mm256_storeu_si256((uint256_t*)(sums_04261537 + 4), sums_1537);

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

				const size_t sum_a = sums_04261537[0];
				*output_a = *input++; // always copy
				// branchless conditional increment
				if ((pf_lookup_ptr[sum_a] & (1 << get_prime_index<17>::idx)) == 0) output_a += 2;

				const size_t sum_b = sums_04261537[4];
				*output_b = *input++;
				if ((pf_lookup_ptr[sum_b] & (1 << get_prime_index<17>::idx)) == 0) output_b += 2;

				const size_t sum_c = sums_04261537[2];
				*output_a = *input++;
				if ((pf_lookup_ptr[sum_c] & (1 << get_prime_index<17>::idx)) == 0) output_a += 2;

				const size_t sum_d = sums_04261537[6];
				*output_b = *input++;
				if ((pf_lookup_ptr[sum_d] & (1 << get_prime_index<17>::idx)) == 0) output_b += 2;

				const size_t sum_e = sums_04261537[1];
				*output_a = *input++;
				if ((pf_lookup_ptr[sum_e] & (1 << get_prime_index<17>::idx)) == 0) output_a += 2;

				const size_t sum_f = sums_04261537[5];
				*output_b = *input++;
				if ((pf_lookup_ptr[sum_f] & (1 << get_prime_index<17>::idx)) == 0) output_b += 2;

				const size_t sum_g = sums_04261537[3];
				*output_a = *input++;
				if ((pf_lookup_ptr[sum_g] & (1 << get_prime_index<17>::idx)) == 0) output_a += 2;

				const size_t sum_h = sums_04261537[7];
				*output_b = *input++;
				if ((pf_lookup_ptr[sum_h] & (1 << get_prime_index<17>::idx)) == 0) output_b += 2;

				// store on the stack for the next iteration
				_mm256_storeu_si256((uint256_t*)sums_04261537, sums_0426);
				_mm256_storeu_si256((uint256_t*)(sums_04261537 + 4), sums_1537);
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
			size_t lookup = prime_factor_lookup[sum] >> get_prime_index<17>::idx;
			output += ~lookup & 0b1;
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
			const prime_lookup_t* const pf_lookup_ptr_12m29 = prime_factor_lookup.data() + b12m29_urem;
			const prime_lookup_t* const pf_lookup_ptr_6m37 = prime_factor_lookup.data() + b6m37_urem;

			uint64_t sums[16]{}; // 0426 1537 0426 1537

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

				// store results on the stack
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_0426_12m29);
				_mm256_storeu_si256((uint256_t*)(sums + 4), sums_1537_12m29);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_0426_6m37);
				_mm256_storeu_si256((uint256_t*)(sums + 12), sums_1537_6m37);
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

				*output = *input++; // always copy
				uint32_t mask_a = 0;
				mask_a |= pf_lookup_ptr_12m29[sums[0]] >> get_prime_index<29>::idx;
				mask_a |= pf_lookup_ptr_6m37[sums[0 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_a & 1; // branchless conditional increment

				*output = *input++;
				uint32_t mask_b = 0;
				mask_b |= pf_lookup_ptr_12m29[sums[4]] >> get_prime_index<29>::idx;
				mask_b |= pf_lookup_ptr_6m37[sums[4 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_b & 1;

				*output = *input++;
				uint32_t mask_c = 0;
				mask_c |= pf_lookup_ptr_12m29[sums[2]] >> get_prime_index<29>::idx;
				mask_c |= pf_lookup_ptr_6m37[sums[2 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_c & 1;

				*output = *input++;
				uint32_t mask_d = 0;
				mask_d |= pf_lookup_ptr_12m29[sums[6]] >> get_prime_index<29>::idx;
				mask_d |= pf_lookup_ptr_6m37[sums[6 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_d & 1;

				*output = *input++;
				uint32_t mask_e = 0;
				mask_e |= pf_lookup_ptr_12m29[sums[1]] >> get_prime_index<29>::idx;
				mask_e |= pf_lookup_ptr_6m37[sums[1 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_e & 1;

				*output = *input++;
				uint32_t mask_f = 0;
				mask_f |= pf_lookup_ptr_12m29[sums[5]] >> get_prime_index<29>::idx;
				mask_f |= pf_lookup_ptr_6m37[sums[5 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_f & 1;

				*output = *input++;
				uint32_t mask_g = 0;
				mask_g |= pf_lookup_ptr_12m29[sums[3]] >> get_prime_index<29>::idx;
				mask_g |= pf_lookup_ptr_6m37[sums[3 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_g & 1;

				*output = *input++;
				uint32_t mask_h = 0;
				mask_h |= pf_lookup_ptr_12m29[sums[7]] >> get_prime_index<29>::idx;
				mask_h |= pf_lookup_ptr_6m37[sums[7 + 8]] >> get_prime_index<37>::idx;
				output += ~mask_h & 1;

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_0426_12m29);
				_mm256_storeu_si256((uint256_t*)(sums + 4), sums_1537_12m29);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_0426_6m37);
				_mm256_storeu_si256((uint256_t*)(sums + 12), sums_1537_6m37);
			} // end for
		} // end if on_fast_path

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
			size_t lookup = 0;
			lookup |= (prime_factor_lookup[sum_12m29] >> get_prime_index<29>::idx);
			lookup |= (prime_factor_lookup[sum_6m37] >> get_prime_index<37>::idx);
			output += ~lookup & 0b1;
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
			const prime_lookup_t* const pf_lookup_ptr_b8m17 = prime_factor_lookup.data() + b8m17_rem;
			const prime_lookup_t* const pf_lookup_ptr_b9m17 = prime_factor_lookup.data() + b9m17_rem;

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
				uint32_t mask_a =
					pf_lookup_ptr_b8m17[sums[0]] |
					pf_lookup_ptr_b9m17[sums[0 + 8]];
				mask_a >>= get_prime_index<17>::idx;
				output += ~mask_a & 1; // branchless conditional increment

				*output = *input++;
				uint32_t mask_b =
					pf_lookup_ptr_b8m17[sums[4]] |
					pf_lookup_ptr_b9m17[sums[4 + 8]];
				mask_b >>= get_prime_index<17>::idx;
				output += ~mask_b & 1;

				*output = *input++;
				uint32_t mask_c =
					pf_lookup_ptr_b8m17[sums[2]] |
					pf_lookup_ptr_b9m17[sums[2 + 8]];
				mask_c >>= get_prime_index<17>::idx;
				output += ~mask_c & 1;

				*output = *input++;
				uint32_t mask_d =
					pf_lookup_ptr_b8m17[sums[6]] |
					pf_lookup_ptr_b9m17[sums[6 + 8]];
				mask_d >>= get_prime_index<17>::idx;
				output += ~mask_d & 1;

				*output = *input++;
				uint32_t mask_e =
					pf_lookup_ptr_b8m17[sums[1]] |
					pf_lookup_ptr_b9m17[sums[1 + 8]];
				mask_e >>= get_prime_index<17>::idx;
				output += ~mask_e & 1;

				*output = *input++;
				uint32_t mask_f =
					pf_lookup_ptr_b8m17[sums[5]] |
					pf_lookup_ptr_b9m17[sums[5 + 8]];
				mask_f >>= get_prime_index<17>::idx;
				output += ~mask_f & 1;

				*output = *input++;
				uint32_t mask_g =
					pf_lookup_ptr_b8m17[sums[3]] |
					pf_lookup_ptr_b9m17[sums[3 + 8]];
				mask_g >>= get_prime_index<17>::idx;
				output += ~mask_g & 1;

				*output = *input++;
				uint32_t mask_h =
					pf_lookup_ptr_b8m17[sums[7]] |
					pf_lookup_ptr_b9m17[sums[7 + 8]];
				mask_h >>= get_prime_index<17>::idx;
				output += ~mask_h & 1;

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_0426_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 4), sums_1537_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_0426_b9m17);
				_mm256_storeu_si256((uint256_t*)(sums + 12), sums_1537_b9m17);
			} // end for
		} // end if on_fast_path

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
			size_t lookup =
				prime_factor_lookup[b8m17_rem] |
				prime_factor_lookup[b9m17_rem];
			lookup >>= get_prime_index<17>::idx;
			output += ~lookup & 0b1;
		}

		return output;

	}

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
