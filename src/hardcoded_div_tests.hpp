#pragma once

#include "config.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/simd.hpp"

namespace mbp
{
	namespace detail
	{
		template<size_t base, size_t prime>
		consteval uint256_t build_4rem_shuffle_lookup()
		{
			uint256_t lookup{ .m256i_u8 = { 0 } };
			for (size_t i = 0; i < 16; ++i)
			{
				if (i & 0x1) lookup.m256i_u8[i] += pow_mod<base, 0, prime>::rem;
				if (i & 0x2) lookup.m256i_u8[i] += pow_mod<base, 1, prime>::rem;
				if (i & 0x4) lookup.m256i_u8[i] += pow_mod<base, 2, prime>::rem;
				if (i & 0x8) lookup.m256i_u8[i] += pow_mod<base, 3, prime>::rem;
				lookup.m256i_u8[i + 16] = lookup.m256i_u8[i]; // duplicate into upper lane
			}
			return lookup;
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

			for (; input < candidates_end_rounded; )
			{
				// load 8 candidates into two registers
				uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));

				// we only want the bottom 32 bits of our 8 candidates
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

				// select the high and low halves of each byte
				const uint256_t nybbles_lo = _mm256_and_si256(ymm0, nybble_mask);
				const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// replace the bits of each nybble with their remainder
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
				merged_masks_0 |= (pf_lookup_ptr_a[sums_0426_a.m256i_u64[0]] >> get_prime_index<5>::idx);
				merged_masks_0 |= (pf_lookup_ptr_b[sums_0426_b.m256i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_0 |= (pf_lookup_ptr_c[sums_0426_c.m256i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_0 |= (pf_lookup_ptr_d[sums_0426_d.m256i_u64[0]] >> get_prime_index<17>::idx);
				*output = *input++; // always copy
				output += ~merged_masks_0 & 0b1; // branchless conditional increment
				size_t merged_masks_4 = 0;
				merged_masks_4 |= (pf_lookup_ptr_a[sums_0426_a.m256i_u64[1]] >> get_prime_index<5>::idx);
				merged_masks_4 |= (pf_lookup_ptr_b[sums_0426_b.m256i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_4 |= (pf_lookup_ptr_c[sums_0426_c.m256i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_4 |= (pf_lookup_ptr_d[sums_0426_d.m256i_u64[1]] >> get_prime_index<17>::idx);

				size_t merged_masks_1 = 0;
				merged_masks_1 |= (pf_lookup_ptr_a[sums_1537_a.m256i_u64[0]] >> get_prime_index<5>::idx);
				merged_masks_1 |= (pf_lookup_ptr_b[sums_1537_b.m256i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_1 |= (pf_lookup_ptr_c[sums_1537_c.m256i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_1 |= (pf_lookup_ptr_d[sums_1537_d.m256i_u64[0]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_1 & 0b1;
				size_t merged_masks_5 = 0;
				merged_masks_5 |= (pf_lookup_ptr_a[sums_1537_a.m256i_u64[1]] >> get_prime_index<5>::idx);
				merged_masks_5 |= (pf_lookup_ptr_b[sums_1537_b.m256i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_5 |= (pf_lookup_ptr_c[sums_1537_c.m256i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_5 |= (pf_lookup_ptr_d[sums_1537_d.m256i_u64[1]] >> get_prime_index<17>::idx);

				uint128_t hi_a = _mm256_extracti128_si256(sums_0426_a, 1);
				uint128_t hi_b = _mm256_extracti128_si256(sums_0426_b, 1);
				uint128_t hi_c = _mm256_extracti128_si256(sums_0426_c, 1);
				uint128_t hi_d = _mm256_extracti128_si256(sums_0426_d, 1);
				size_t merged_masks_2 = 0;
				merged_masks_2 |= (pf_lookup_ptr_a[hi_a.m128i_u64[0]] >> get_prime_index<5>::idx);
				merged_masks_2 |= (pf_lookup_ptr_b[hi_b.m128i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_2 |= (pf_lookup_ptr_c[hi_c.m128i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_2 |= (pf_lookup_ptr_d[hi_d.m128i_u64[0]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_2 & 0b1;
				size_t merged_masks_6 = 0;
				merged_masks_6 |= (pf_lookup_ptr_a[hi_a.m128i_u64[1]] >> get_prime_index<5>::idx);
				merged_masks_6 |= (pf_lookup_ptr_b[hi_b.m128i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_6 |= (pf_lookup_ptr_c[hi_c.m128i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_6 |= (pf_lookup_ptr_d[hi_d.m128i_u64[1]] >> get_prime_index<17>::idx);

				hi_a = _mm256_extracti128_si256(sums_1537_a, 1);
				hi_b = _mm256_extracti128_si256(sums_1537_b, 1);
				hi_c = _mm256_extracti128_si256(sums_1537_c, 1);
				hi_d = _mm256_extracti128_si256(sums_1537_d, 1);
				size_t merged_masks_3 = 0;
				merged_masks_3 |= (pf_lookup_ptr_a[hi_a.m128i_u64[0]] >> get_prime_index<5>::idx);
				merged_masks_3 |= (pf_lookup_ptr_b[hi_b.m128i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_3 |= (pf_lookup_ptr_c[hi_c.m128i_u64[0]] >> get_prime_index<13>::idx);
				merged_masks_3 |= (pf_lookup_ptr_d[hi_d.m128i_u64[0]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_3 & 0b1;
				*output = *input++;
				output += ~merged_masks_4 & 0b1;
				*output = *input++;
				output += ~merged_masks_5 & 0b1;
				*output = *input++;
				output += ~merged_masks_6 & 0b1;

				size_t merged_masks_7 = 0;
				merged_masks_7 |= (pf_lookup_ptr_a[hi_a.m128i_u64[1]] >> get_prime_index<5>::idx);
				merged_masks_7 |= (pf_lookup_ptr_b[hi_b.m128i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_7 |= (pf_lookup_ptr_c[hi_c.m128i_u64[1]] >> get_prime_index<13>::idx);
				merged_masks_7 |= (pf_lookup_ptr_d[hi_d.m128i_u64[1]] >> get_prime_index<17>::idx);
				*output = *input++;
				output += ~merged_masks_7 & 0b1;
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

		// constexpr uint128_t static_rems0 = { always all ones };
		constexpr static uint128_t static_rems1 = uint128_t{ .m128i_u32{
			pow_mod<3, 1, 7>::rem,
			pow_mod<5, 1, 7>::rem,
			pow_mod<4, 1, 13>::rem,
			pow_mod<10, 1, 13>::rem } };
		constexpr static uint128_t static_rems2 = uint128_t{ .m128i_u32{
			pow_mod<3, 2, 7>::rem,
			pow_mod<5, 2, 7>::rem,
			pow_mod<4, 2, 13>::rem,
			pow_mod<10, 2, 13>::rem } };
		constexpr static uint128_t static_rems3 = uint128_t{ .m128i_u32{
			pow_mod<3, 3, 7>::rem,
			pow_mod<5, 3, 7>::rem,
			pow_mod<4, 3, 13>::rem,
			pow_mod<10, 3, 13>::rem } };
		constexpr static uint128_t static_rems4 = uint128_t{ .m128i_u32{
			pow_mod<3, 4, 7>::rem,
			pow_mod<5, 4, 7>::rem,
			pow_mod<4, 4, 13>::rem,
			pow_mod<10, 4, 13>::rem } };
		constexpr static uint128_t static_rems5 = uint128_t{ .m128i_u32{
			pow_mod<3, 5, 7>::rem,
			pow_mod<5, 5, 7>::rem,
			pow_mod<4, 5, 13>::rem,
			pow_mod<10, 5, 13>::rem } };

		const auto xmm_rems1 = _mm_loadu_si128(&static_rems1);
		const auto xmm_rems2 = _mm_loadu_si128(&static_rems2);
		const auto xmm_rems3 = _mm_loadu_si128(&static_rems3);
		const auto xmm_rems4 = _mm_loadu_si128(&static_rems4);
		const auto xmm_rems5 = _mm_loadu_si128(&static_rems5);

		const prime_lookup_t* const prime_factor_lookup_ptr = prime_factor_lookup.data();

		size_t* output = input;

		size_t next = *input;

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			auto xmm0 = _mm_set1_epi32((int)pop_count(number & (bitmask << 0)));
			auto xmm1 = _mm_set1_epi32((int)pop_count(number & (bitmask << 1)));
			auto xmm2 = _mm_set1_epi32((int)pop_count(number & (bitmask << 2)));
			auto xmm3 = _mm_set1_epi32((int)pop_count(number & (bitmask << 3)));
			auto xmm4 = _mm_set1_epi32((int)pop_count(number & (bitmask << 4)));
			auto xmm5 = _mm_set1_epi32((int)pop_count(number & (bitmask << 5)));

			// multiply each popcount with Nth remainder of each test
			// (rems0 is always all ones)
			xmm1 = _mm_mullo_epi32(xmm1, xmm_rems1);
			xmm2 = _mm_mullo_epi32(xmm2, xmm_rems2);
			xmm3 = _mm_mullo_epi32(xmm3, xmm_rems3);
			xmm4 = _mm_mullo_epi32(xmm4, xmm_rems4);
			xmm5 = _mm_mullo_epi32(xmm5, xmm_rems5);

			// add pairs
			xmm0 = _mm_add_epi32(xmm0, xmm1);
			xmm2 = _mm_add_epi32(xmm2, xmm3);
			xmm4 = _mm_add_epi32(xmm4, xmm5);
			// final adds
			xmm0 = _mm_add_epi32(xmm0, xmm2);
			xmm0 = _mm_add_epi32(xmm0, xmm4);

			// Only advance the pointer if the number is still a candidate
			size_t merged_masks = 0;
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[0]] >> get_prime_index<7>::idx);
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[1]] >> get_prime_index<7>::idx);
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[2]] >> get_prime_index<13>::idx);
			merged_masks |= (prime_factor_lookup_ptr[xmm0.m128i_u32[3]] >> get_prime_index<13>::idx);

			output += ~merged_masks & 0b1;
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

		const int* const prime_factor_lookup_ptr = (const int*)prime_factor_lookup.data();

		size_t* output = input;

		for (; input < candidates_end; ++input)
		{
			const size_t number = *input;

			// always write
			*output = number;

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

			// Only advance the pointer if the number is still a candidate
			size_t merged_lookups = prime_factor_lookup_ptr[xmm0.m128i_u32[0]];
			merged_lookups |= prime_factor_lookup_ptr[xmm0.m128i_u32[1]];
			merged_lookups |= prime_factor_lookup_ptr[xmm0.m128i_u32[2]];
			merged_lookups |= prime_factor_lookup_ptr[xmm0.m128i_u32[3]];

			merged_lookups >>= get_prime_index<11>::idx;

			output += ~merged_lookups & 0b1;
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

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		static_assert(sizeof(remainder_t) == 1);
		alignas(64) static constexpr remainder_t b6_rems[64] = {
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7, 9, 10, 5, 8, 4, 2,
			1, 6, 3, 7 };
		alignas(64) static constexpr remainder_t b7_rems[64] = {
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2, 3, 10, 4, 6, 9, 8,
			1, 7, 5, 2 };
		alignas(64) static constexpr remainder_t b8_rems[64] = {
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6, 4, 10, 3, 2, 5, 7,
			1, 8, 9, 6 };

		const uint256_t b6_rems_lower = _mm256_loadu_si256((uint256_t*)&b6_rems[0]);
		const uint256_t b6_rems_upper = _mm256_loadu_si256((uint256_t*)&b6_rems[32]);

		const uint256_t b7_rems_lower = _mm256_loadu_si256((uint256_t*)&b7_rems[0]);
		const uint256_t b7_rems_upper = _mm256_loadu_si256((uint256_t*)&b7_rems[32]);

		const uint256_t b8_rems_lower = _mm256_loadu_si256((uint256_t*)&b8_rems[0]);
		const uint256_t b8_rems_upper = _mm256_loadu_si256((uint256_t*)&b8_rems[32]);

		size_t* output = input;
		size_t next = *input;

		size_t upper_bits = next & upper_bits_mask;
		uint256_t mask_upper = util::expand_bits_to_bytes(next >> 32);
		uint256_t ymm1 = _mm256_and_si256(mask_upper, b6_rems_upper);
		uint256_t ymm3 = _mm256_and_si256(mask_upper, b7_rems_upper);
		uint256_t ymm5 = _mm256_and_si256(mask_upper, b8_rems_upper);

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			if constexpr (!on_fast_path)
			{
				// Upper 32 bits only change every 4B integers - re-calculate when stale
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;
					mask_upper = util::expand_bits_to_bytes(number >> 32);
					ymm1 = _mm256_and_si256(mask_upper, b6_rems_upper);
					ymm3 = _mm256_and_si256(mask_upper, b7_rems_upper);
					ymm5 = _mm256_and_si256(mask_upper, b8_rems_upper);
				}
			}

			// Convert 64 bits to 64 bytes
			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));

			const uint256_t ymm0 = _mm256_and_si256(mask_lower, b6_rems_lower);
			const uint256_t ymm2 = _mm256_and_si256(mask_lower, b7_rems_lower);
			const uint256_t ymm4 = _mm256_and_si256(mask_lower, b8_rems_lower);

			// Not all div tests can safely combine rems without moving to a larger type,
			// but our rems are small enough to avoid overflow.
			const auto b6_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm0, ymm1));
			const auto b7_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm2, ymm3));
			const auto b8_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm4, ymm5));

			// Only advance the pointer if the number is still a candidate
			size_t merged_lookups = prime_factor_lookup[b6_rem];
			merged_lookups |= prime_factor_lookup[b7_rem];
			merged_lookups |= prime_factor_lookup[b8_rem];

			merged_lookups >>= get_prime_index<11>::idx;

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

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		static_assert(sizeof(remainder_t) == 1);
		alignas(64) static constexpr remainder_t b6_rems[64] = {
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11,
			1, 6, 10, 8 };
		alignas(64) static constexpr remainder_t b7_rems[64] = {
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2,
			1, 7, 10, 5 };
		alignas(64) static constexpr remainder_t b11_rems[64] = {
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6,
			1, 11, 4, 5 };

		const uint256_t b6_rems_lower = _mm256_loadu_si256((uint256_t*)&b6_rems[0]);
		const uint256_t b6_rems_upper = _mm256_loadu_si256((uint256_t*)&b6_rems[32]);

		const uint256_t b7_rems_lower = _mm256_loadu_si256((uint256_t*)&b7_rems[0]);
		const uint256_t b7_rems_upper = _mm256_loadu_si256((uint256_t*)&b7_rems[32]);

		const uint256_t b11_rems_lower = _mm256_loadu_si256((uint256_t*)&b11_rems[0]);
		const uint256_t b11_rems_upper = _mm256_loadu_si256((uint256_t*)&b11_rems[32]);

		size_t* output = input;
		size_t next = *input;

		size_t upper_bits = next & upper_bits_mask;
		uint256_t mask_upper = util::expand_bits_to_bytes(next >> 32);
		uint256_t ymm1 = _mm256_and_si256(mask_upper, b6_rems_upper);
		uint256_t ymm3 = _mm256_and_si256(mask_upper, b7_rems_upper);
		uint256_t ymm5 = _mm256_and_si256(mask_upper, b11_rems_upper);

		for (; input < candidates_end; )
		{
			const size_t number = next;
			++input;
			next = *input; // load one iteration ahead

			// always write
			*output = number;

			if constexpr (!on_fast_path)
			{
				// Upper 32 bits only change every 4B integers - re-calculate when stale
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;
					mask_upper = util::expand_bits_to_bytes(number >> 32);
					ymm1 = _mm256_and_si256(mask_upper, b6_rems_upper);
					ymm3 = _mm256_and_si256(mask_upper, b7_rems_upper);
					ymm5 = _mm256_and_si256(mask_upper, b11_rems_upper);
				}
			}

			const uint256_t mask_lower = util::expand_bits_to_bytes(number & uint32_t(-1));

			const uint256_t ymm0 = _mm256_and_si256(mask_lower, b6_rems_lower);
			const uint256_t ymm2 = _mm256_and_si256(mask_lower, b7_rems_lower);
			const uint256_t ymm4 = _mm256_and_si256(mask_lower, b11_rems_lower);

			// Not all div tests can safely combine rems without moving to a larger type,
			// but our rems are small enough to avoid overflow.
			const auto b6_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm0, ymm1));
			const auto b7_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm2, ymm3));
			const auto b11_rem = util::vcl_hadd_x(_mm256_add_epi8(ymm4, ymm5));

			// Only advance the pointer if the number is still a candidate
			prime_lookup_t merged_lookups = prime_factor_lookup[b6_rem];
			merged_lookups |= prime_factor_lookup[b7_rem];
			merged_lookups |= prime_factor_lookup[b11_rem];

			merged_lookups >>= get_prime_index<13>::idx;
			
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

		constexpr static uint8_t alignas(64) static_rems[7 + 1][16] = {
			{ 1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6 }, // base 3 % 17
			{ 1, 5, 8, 6, 13, 14, 2, 10, 16, 12, 9, 11, 4, 3, 15, 7 }, // base 5 % 17
			{ 1, 6, 2, 12, 4, 7, 8, 14, 16, 11, 15, 5, 13, 10, 9, 3 }, // base 6 % 17
			{ 1, 7, 15, 3, 4, 11, 9, 12, 16, 10, 2, 14, 13, 6, 8, 5 }, // base 7 % 17
			{ 1, 10, 15, 14, 4, 6, 9, 5, 16, 7, 2, 3, 13, 11, 8, 12 }, // base 10 % 17
			{ 1, 11, 2, 5, 4, 10, 8, 3, 16, 6, 15, 12, 13, 7, 9, 14 }, // base 11 % 17
			{ 1, 12, 8, 11, 13, 3, 2, 7, 16, 5, 9, 6, 4, 14, 15, 10 }, // base 12 % 17
			{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } }; // padding

		const uint256_t rems_0_1 = _mm256_loadu_si256((uint256_t*)&static_rems[0]);
		const uint256_t rems_2_3 = _mm256_loadu_si256((uint256_t*)&static_rems[2]);
		const uint256_t rems_4_5 = _mm256_loadu_si256((uint256_t*)&static_rems[4]);
		const uint256_t rems_6_x = _mm256_loadu_si256((uint256_t*)&static_rems[6]);

		size_t* output = input;
		size_t number = *input;

		size_t upper_bits = number & upper_bits_mask;
		uint256_t upper_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number >> 32),
												 _mm256_set1_epi8(0x01));

		for (; input < candidates_end; )
		{
			if constexpr (!on_fast_path)
			{
				// Upper 32 bits only change every 4B integers - re-calculate when stale
				if ((number & upper_bits_mask) != upper_bits)
				{
					upper_bits = number & upper_bits_mask;
					upper_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number >> 32),
												   _mm256_set1_epi8(0x01)); // convert 0xFF -> 0x01
				}
			}

			uint256_t lower_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number & uint32_t(-1)),
													 _mm256_set1_epi8(0x01)); // convert 0xFF -> 0x01

			// always write
			*output = number;
			// load one iteration ahead
			++input;
			number = *input;

			// add the upper 32 bytes to the lower 32 bytes to get values in range 0-2
			uint256_t ymm0 = _mm256_add_epi8(lower_bytes, upper_bytes);

			// add the upper 16 bytes to the lower 16 bytes to get values in range 0-4, and leave the results in both lanes
			ymm0 = _mm256_add_epi8(ymm0, _mm256_permute2x128_si256(ymm0, ymm0, 1));

			// multiply 16+16 8-bit values by 16+16 remainders, storing partially summed results as 8+8 uint16_ts
			uint256_t ymm1 = _mm256_maddubs_epi16(ymm0, rems_0_1);
			uint256_t ymm2 = _mm256_maddubs_epi16(ymm0, rems_2_3);
			uint256_t ymm3 = _mm256_maddubs_epi16(ymm0, rems_4_5);
			uint256_t ymm4 = _mm256_maddubs_epi16(ymm0, rems_6_x);

			// pack four sparse registers into two
			ymm0 = _mm256_packus_epi16(ymm1, ymm2);
			ymm1 = _mm256_packus_epi16(ymm3, ymm4);

			// h-sum into 4 16-bit integers, 2 in each lane
			ymm0 = _mm256_sad_epu8(ymm0, _mm256_setzero_si256());
			ymm1 = _mm256_sad_epu8(ymm1, _mm256_setzero_si256());

			// extract high lanes
			const uint128_t xmm2 = _mm256_extracti128_si256(ymm0, 1);
			const uint128_t xmm3 = _mm256_extracti128_si256(ymm1, 1);

			prime_lookup_t merged_lookups = 0;
			merged_lookups |= prime_factor_lookup[ymm0.m256i_i64[0]];
			merged_lookups |= prime_factor_lookup[ymm0.m256i_i64[1]];
			merged_lookups |= prime_factor_lookup[ymm1.m256i_i64[0]];
			merged_lookups |= prime_factor_lookup[ymm1.m256i_i64[1]];
			merged_lookups |= prime_factor_lookup[xmm2.m128i_i64[0]];
			merged_lookups |= prime_factor_lookup[xmm2.m128i_i64[1]];
			merged_lookups |= prime_factor_lookup[xmm3.m128i_i64[0]];

			// Only advance the pointer if the nth bit was 0 in all lookups
			// Invert the bit we care about, move it to the 0th position, and mask to select it
			output += (~merged_lookups >> get_prime_index<17>::idx) & 0b1;
		}

		return output;
	}

	// currently unused:
	template<bool on_fast_path>
	tests_are_inlined size_t* div_tests_with_8_rems(size_t* input,
													const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		// base 12 % 29: 2*4 remainders: 1  12  28  17
		// base  9 % 17:   8 remainders: 1   9  13  15   16   8   4   2
		// base  6 % 37: 2*4 remainders: 1   6  36  31
		// base  8 % 17:   8 remainders: 1   8  13   2   16   9   4  15
		constexpr static uint256_t static_rems_0_1 = { .m256i_u8{
			1, 12, 28, 17, 1, 12, 28, 17,
			1, 9, 13, 15, 16, 8, 4, 2,
			1, 12, 28, 17, 1, 12, 28, 17,
			1, 9, 13, 15, 16, 8, 4, 2 } };
		constexpr static uint256_t static_rems_2_3 = { .m256i_u8{
			1, 6, 36, 31, 1, 6, 36, 31,
			1, 8, 13, 2, 16, 9, 4, 15,
			1, 6, 36, 31, 1, 6, 36, 31,
			1, 8, 13, 2, 16, 9, 4, 15 } };

		constexpr static uint256_t static_bits_to_bytes_shuffle_mask = { .m256i_u8{
			0, 0, 1, 1, 2, 2, 3, 3,
			0, 0, 1, 1, 2, 2, 3, 3,
			8, 8, 9, 9, 10, 10, 11, 11, // access eight bytes higher in the upper lane
			8, 8, 9, 9, 10, 10, 11, 11 } };
		constexpr static uint256_t static_and_mask = { .m256i_u64 = {
			0x10'01'10'01'10'01'10'01,
			0x20'02'20'02'20'02'20'02,
			0x10'01'10'01'10'01'10'01,
			0x20'02'20'02'20'02'20'02 } };

		const uint256_t rems_0_1 = _mm256_loadu_si256(&static_rems_0_1);
		const uint256_t rems_2_3 = _mm256_loadu_si256(&static_rems_2_3);

		const uint256_t bits_to_bytes_shuffle_mask = _mm256_loadu_si256(&static_bits_to_bytes_shuffle_mask);
		const uint256_t and_mask = _mm256_loadu_si256(&static_and_mask);
		const uint256_t one_bit_mask = _mm256_set1_epi8(0x1);

		size_t* output = input;

		size_t number_a = *input;
		size_t number_b = *(input + 1);

		size_t upper_bits_a = number_a & upper_bits_mask;
		size_t upper_bits_b = number_b & upper_bits_mask;

		uint256_t upper_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number_a >> 32),
												 one_bit_mask);

		for (; input < candidates_end; input += 2)
		{
			const size_t number_a = *input;
			const size_t number_b = *(input + 1);

			if constexpr (!on_fast_path)
			{
				// todo

				// Upper 32 bits only change every 4B integers - re-calculate when stale
				if ((number_a & upper_bits_mask) != upper_bits_a)
				{
					upper_bits_a = number_a & upper_bits_mask;
					upper_bytes = _mm256_and_si256(util::expand_bits_to_bytes(number_a >> 32),
												   one_bit_mask); // convert 0xFF -> 0x01
				}
			}

			// load two candidates
			const uint128_t xmm_numbers = _mm_loadu_si128((uint128_t*)input);
			const uint256_t numbers = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_numbers), xmm_numbers, 1);

			// Convert the low 16 bits of each candidate to 16 bytes, and repeat for the next 16 bits
			uint256_t lo_two_bytes = _mm256_shuffle_epi8(numbers, bits_to_bytes_shuffle_mask);
			uint256_t hi_two_bytes = _mm256_srli_epi64(numbers, 2); // shift by 2 bits
			hi_two_bytes = _mm256_shuffle_epi8(hi_two_bytes, bits_to_bytes_shuffle_mask);
			lo_two_bytes = _mm256_andnot_si256(lo_two_bytes, and_mask);
			hi_two_bytes = _mm256_andnot_si256(hi_two_bytes, and_mask);
			lo_two_bytes = _mm256_cmpeq_epi8(lo_two_bytes, _mm256_setzero_si256());
			hi_two_bytes = _mm256_cmpeq_epi8(hi_two_bytes, _mm256_setzero_si256());
			lo_two_bytes = _mm256_and_si256(lo_two_bytes, one_bit_mask);
			hi_two_bytes = _mm256_and_si256(hi_two_bytes, one_bit_mask);

			// add bytes to get 16+16 values in the range of 0-2
			uint256_t lower_bytes = _mm256_add_epi8(lo_two_bytes, hi_two_bytes);
			// swap upper and lower half of each lane, and add to get values in range 0-4
			lower_bytes = _mm256_add_epi8(lower_bytes, _mm256_shuffle_epi32(lower_bytes, 0b01'00'11'10));

			// we could add lower 4 and upper 4 bytes here, but I think that can be done later

			// multiply 2 sets of 2 copies of 8 rems with 2 tests, then 2 more tests
			uint256_t sums_0_1 = _mm256_maddubs_epi16(lower_bytes, rems_0_1);
			uint256_t sums_2_3 = _mm256_maddubs_epi16(lower_bytes, rems_2_3);

			// h-sum partial sums to 4x sums
			sums_0_1 = _mm256_sad_epu8(sums_0_1, _mm256_setzero_si256());
			sums_2_3 = _mm256_sad_epu8(sums_2_3, _mm256_setzero_si256());




			// add the upper 32 bytes to the lower 32 bytes to get values in range 0-2
			uint256_t ymm0 = _mm256_add_epi8(lower_bytes, upper_bytes);

			// multiply 16+16 8-bit values by 32 remainders, storing partially summed results as 8+8 uint16_ts
			uint256_t ymm1 = _mm256_maddubs_epi16(ymm0, rems_0_1);

			// [ some final add step(s) ]

			// extract high lanes
			const uint128_t xmm1 = _mm256_extracti128_si256(ymm0, 1);

			prime_lookup_t merged_lookups_a = 0;
			merged_lookups_a |= prime_factor_lookup[ymm0.m256i_i64[0]];
			merged_lookups_a |= prime_factor_lookup[ymm0.m256i_i64[1]];
			merged_lookups_a |= prime_factor_lookup[xmm1.m128i_i64[0]];
			merged_lookups_a |= prime_factor_lookup[xmm1.m128i_i64[1]];

			prime_lookup_t merged_lookups_b = 0;

			// Only advance the pointer if the number is still a candidate, that is,
			// the relevant bit from each lookup is 0

			*output = number_a;
			output += ((~merged_lookups_a) & 0b1);

			*output = number_b;
			output += ((~merged_lookups_b) & 0b1);
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
				uint32_t mask_a = pf_lookup_ptr[sums_0426.m256i_u64[0]] >> get_prime_index<17>::idx;
				output += ~mask_a & 1; // branchless conditional increment

				*output = *input++;
				uint32_t mask_b = pf_lookup_ptr[sums_1537.m256i_u64[0]] >> get_prime_index<17>::idx;
				output += ~mask_b & 1;

				*output = *input++;
				uint32_t mask_c = pf_lookup_ptr[sums_26.m128i_u64[0]] >> get_prime_index<17>::idx;
				output += ~mask_c & 1;

				*output = *input++;
				uint32_t mask_d = pf_lookup_ptr[sums_37.m128i_u64[0]] >> get_prime_index<17>::idx;
				output += ~mask_d & 1;

				*output = *input++;
				uint32_t mask_e = pf_lookup_ptr[sums_0426.m256i_u64[1]] >> get_prime_index<17>::idx;
				output += ~mask_e & 1;

				*output = *input++;
				uint32_t mask_f = pf_lookup_ptr[sums_1537.m256i_u64[1]] >> get_prime_index<17>::idx;
				output += ~mask_f & 1;

				*output = *input++;
				uint32_t mask_g = pf_lookup_ptr[sums_26.m128i_u64[1]] >> get_prime_index<17>::idx;
				output += ~mask_g & 1;

				*output = *input++;
				uint32_t mask_h = pf_lookup_ptr[sums_37.m128i_u64[1]] >> get_prime_index<17>::idx;
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
}
