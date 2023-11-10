#pragma once

#include "config.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/simd.hpp"

namespace mbp
{
	namespace detail
	{
		template<size_t base, size_t place_value_start, size_t prime>
		consteval std::array<uint8_t, 32> build_shuffle_lookup_impl()
		{
			std::array<uint8_t, 32> lookup{};
			for (size_t i = 0; i < 16; ++i)
			{
				if (i & 0b0001) lookup[i] += pow_mod(base, place_value_start + 0, prime);
				if (i & 0b0010) lookup[i] += pow_mod(base, place_value_start + 1, prime);
				if (i & 0b0100) lookup[i] += pow_mod(base, place_value_start + 2, prime);
				if (i & 0b1000) lookup[i] += pow_mod(base, place_value_start + 3, prime);
				lookup[i + 16] = lookup[i]; // duplicate into upper lane
			}
			return lookup;
		}

		template<size_t base, size_t prime>
		consteval std::array<uint8_t, 32> build_4rem_shuffle_lookup()
		{
			return build_shuffle_lookup_impl<base, 0, prime>();
		}

		template<size_t base, size_t prime>
		consteval std::array<uint8_t, 32> build_8rem_shuffle_lookup_lo_nybble()
		{
			return build_shuffle_lookup_impl<base, 0, prime>();
		}

		template<size_t base, size_t prime>
		consteval std::array<uint8_t, 32> build_8rem_shuffle_lookup_hi_nybble()
		{
			return build_shuffle_lookup_impl<base, 4, prime>();
		}
	}

	template<size_t base_a, size_t prime_a, size_t base_b, size_t prime_b, bool on_fast_path>
	inline_toggle size_t* two_div_tests_with_four_rems(size_t* input,
														   const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t bitmask = bitmask_for<base_a, prime_a>::val;
		static_assert(bitmask == bitmask_for<base_b, prime_b>::val);
		static_assert(period_of<bitmask>() == 4);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static std::array<uint8_t, 32> static_nybble_lookup_a = mbp::detail::build_4rem_shuffle_lookup<base_a, prime_a>();
			constexpr static std::array<uint8_t, 32> static_nybble_lookup_b = mbp::detail::build_4rem_shuffle_lookup<base_b, prime_b>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup_a = _mm256_loadu_si256((uint256_t*)&static_nybble_lookup_a);
			const uint256_t nybble_lookup_b = _mm256_loadu_si256((uint256_t*)&static_nybble_lookup_b);

			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			const size_t number_upper = (*input) & (size_t(-1) << 32);
			const size_t pc_0 = pop_count(number_upper & (bitmask << 0));
			size_t sum_a = pc_0;
			size_t sum_b = pc_0;
			const size_t pc_1 = pop_count(number_upper & (bitmask << 1));
			sum_a += pc_1 * pow_mod(base_a, 1, prime_a);
			sum_b += pc_1 * pow_mod(base_b, 1, prime_b);
			const size_t pc_2 = pop_count(number_upper & (bitmask << 2));
			sum_a += pc_2 * pow_mod(base_a, 2, prime_a);
			sum_b += pc_2 * pow_mod(base_b, 2, prime_b);
			const size_t pc_3 = pop_count(number_upper & (bitmask << 3));
			sum_a += pc_3 * pow_mod(base_a, 3, prime_a);
			sum_b += pc_3 * pow_mod(base_b, 3, prime_b);
			const uint8_t* const indivisible_base_a = indivisible_by[get_prime_index<prime_a>::idx].data() + sum_a;
			const uint8_t* const indivisible_base_b = indivisible_by[get_prime_index<prime_b>::idx].data() + sum_b;

			alignas(32) volatile uint32_t sums[16]{};

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
				const uint256_t lo_nybbles_a = _mm256_shuffle_epi8(nybble_lookup_a, lo_nybbles);
				const uint256_t hi_nybbles_a = _mm256_shuffle_epi8(nybble_lookup_a, hi_nybbles);
				const uint256_t lo_nybbles_b = _mm256_shuffle_epi8(nybble_lookup_b, lo_nybbles);
				const uint256_t hi_nybbles_b = _mm256_shuffle_epi8(nybble_lookup_b, hi_nybbles);
				// repack
				uint256_t sums_0426_a = _mm256_unpacklo_epi32(lo_nybbles_a, hi_nybbles_a);
				uint256_t sums_1537_a = _mm256_unpackhi_epi32(lo_nybbles_a, hi_nybbles_a);
				uint256_t sums_0426_b = _mm256_unpacklo_epi32(lo_nybbles_b, hi_nybbles_b);
				uint256_t sums_1537_b = _mm256_unpackhi_epi32(lo_nybbles_b, hi_nybbles_b);
				// h-sum each set of 8 remainders
				sums_0426_a = _mm256_sad_epu8(sums_0426_a, _mm256_setzero_si256());
				sums_1537_a = _mm256_sad_epu8(sums_1537_a, _mm256_setzero_si256());
				sums_0426_b = _mm256_sad_epu8(sums_0426_b, _mm256_setzero_si256());
				sums_1537_b = _mm256_sad_epu8(sums_1537_b, _mm256_setzero_si256());

				const uint256_t sums_0404_2626 = _mm256_packus_epi32(sums_0426_a, sums_0426_b);
				const uint256_t sums_1515_3737 = _mm256_packus_epi32(sums_1537_a, sums_1537_b);

				// store results on the stack
				_mm256_storeu_si256((uint256_t*)&sums[0], sums_0404_2626);
				_mm256_storeu_si256((uint256_t*)&sums[8], sums_1515_3737);
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

				const uint256_t lo_nybbles_a = _mm256_shuffle_epi8(nybble_lookup_a, lo_nybbles);
				const uint256_t hi_nybbles_a = _mm256_shuffle_epi8(nybble_lookup_a, hi_nybbles);
				const uint256_t lo_nybbles_b = _mm256_shuffle_epi8(nybble_lookup_b, lo_nybbles);
				const uint256_t hi_nybbles_b = _mm256_shuffle_epi8(nybble_lookup_b, hi_nybbles);
				uint256_t sums_0426_a = _mm256_unpacklo_epi32(lo_nybbles_a, hi_nybbles_a);
				uint256_t sums_1537_a = _mm256_unpackhi_epi32(lo_nybbles_a, hi_nybbles_a);
				uint256_t sums_0426_b = _mm256_unpacklo_epi32(lo_nybbles_b, hi_nybbles_b);
				uint256_t sums_1537_b = _mm256_unpackhi_epi32(lo_nybbles_b, hi_nybbles_b);
				sums_0426_a = _mm256_sad_epu8(sums_0426_a, _mm256_setzero_si256());
				sums_1537_a = _mm256_sad_epu8(sums_1537_a, _mm256_setzero_si256());
				sums_0426_b = _mm256_sad_epu8(sums_0426_b, _mm256_setzero_si256());
				sums_1537_b = _mm256_sad_epu8(sums_1537_b, _mm256_setzero_si256());

				const uint256_t sums_0404_2626 = _mm256_packus_epi32(sums_0426_a, sums_0426_b);
				const uint256_t sums_1515_3737 = _mm256_packus_epi32(sums_1537_a, sums_1537_b);

				*output = *input++; // always copy
				const size_t inc_a =
					indivisible_base_a[sums[0]] &
					indivisible_base_b[sums[0 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_a); // branchless conditional increment

				*output = *input++;
				const size_t inc_b =
					indivisible_base_a[sums[8]] &
					indivisible_base_b[sums[8 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_b);

				*output = *input++;
				const size_t inc_c =
					indivisible_base_a[sums[4]] &
					indivisible_base_b[sums[4 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_c);

				*output = *input++;
				const size_t inc_d =
					indivisible_base_a[sums[12]] &
					indivisible_base_b[sums[12 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_d);

				*output = *input++;
				const size_t inc_e =
					indivisible_base_a[sums[1]] &
					indivisible_base_b[sums[1 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_e);

				*output = *input++;
				const size_t inc_f =
					indivisible_base_a[sums[9]] &
					indivisible_base_b[sums[9 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_f);

				*output = *input++;
				const size_t inc_g =
					indivisible_base_a[sums[5]] &
					indivisible_base_b[sums[5 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_g);

				*output = *input++;
				const size_t inc_h =
					indivisible_base_a[sums[13]] &
					indivisible_base_b[sums[13 + 2]];
				output = (uint64_t*)(((uint8_t*)output) + inc_h);

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)&sums[0], sums_0404_2626);
				_mm256_storeu_si256((uint256_t*)&sums[8], sums_1515_3737);
			} // end for
		} // end if on_fast_path

		// handle any elements not handled by the fast loop

		const uint8_t* const indivisible_base_a = indivisible_by[get_prime_index<prime_a>::idx].data();
		const uint8_t* const indivisible_base_b = indivisible_by[get_prime_index<prime_b>::idx].data();

		size_t candidate = *input;
		for (; input < candidates_end; )
		{
			*output = candidate; // always write

			const size_t pc_0 = pop_count(candidate & (bitmask << 0));
			size_t sum_a = pc_0;
			size_t sum_b = pc_0;
			const size_t pc_1 = pop_count(candidate & (bitmask << 1));
			sum_a += pc_1 * pow_mod(base_a, 1, prime_a);
			sum_b += pc_1 * pow_mod(base_b, 1, prime_b);
			const size_t pc_2 = pop_count(candidate & (bitmask << 2));
			sum_a += pc_2 * pow_mod(base_a, 2, prime_a);
			sum_b += pc_2 * pow_mod(base_b, 2, prime_b);
			const size_t pc_3 = pop_count(candidate & (bitmask << 3));
			sum_a += pc_3 * pow_mod(base_a, 3, prime_a);
			sum_b += pc_3 * pow_mod(base_b, 3, prime_b);

			// load one iteration ahead
			candidate = *++input;

			// conditionally increment
			const size_t inc =
				indivisible_base_a[sum_a] &
				indivisible_base_b[sum_b];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
	inline_toggle size_t* div_tests_with_12_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		// bases 6, 7, and 11 mod 13 (12 remainders)
		constexpr size_t bitmask = bitmask_for<6, 13>::val;
		static_assert(bitmask == bitmask_for<7, 13>::val);
		static_assert(bitmask == bitmask_for<11, 13>::val);
		static_assert(period_of<bitmask>() == 12);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr size_t upper_bits_mask = size_t(-1) << 32;

			constexpr static std::array<uint8_t, 16> static_b6_rems_lo = { 1, 6, 10, 8, 9, 2, 12, 7, 3, 5, 4, 11, 1, 6, 10, 8 };
			constexpr static std::array<uint8_t, 16> static_b6_rems_hi = { 9, 2, 12, 7, 3, 5, 4, 11, 1, 6, 10, 8, 9, 2, 12, 7 };
			constexpr static std::array<uint8_t, 16> static_b7_rems_lo = { 1, 7, 10, 5, 9, 11, 12, 6, 3, 8, 4, 2, 1, 7, 10, 5 };
			constexpr static std::array<uint8_t, 16> static_b7_rems_hi = { 9, 11, 12, 6, 3, 8, 4, 2, 1, 7, 10, 5, 9, 11, 12, 6 };
			constexpr static std::array<uint8_t, 16> static_b11_rems_lo = { 1, 11, 4, 5, 3, 7, 12, 2, 9, 8, 10, 6, 1, 11, 4, 5 };
			constexpr static std::array<uint8_t, 16> static_b11_rems_hi = { 3, 7, 12, 2, 9, 8, 10, 6, 1, 11, 4, 5, 3, 7, 12, 2 };

			const size_t upper_bits = (*input) & upper_bits_mask;
			const size_t pc_0 = pop_count(upper_bits & (bitmask << 0));
			size_t b6_sum = pc_0;
			size_t b7_sum = pc_0;
			size_t b11_sum = pc_0;
			const size_t pc_1 = pop_count(upper_bits & (bitmask << 1));
			b6_sum += pc_1 * pow_mod(6, 1, 13);
			b7_sum += pc_1 * pow_mod(7, 1, 13);
			b11_sum += pc_1 * pow_mod(11, 1, 13);
			const size_t pc_2 = pop_count(upper_bits & (bitmask << 2));
			b6_sum += pc_2 * pow_mod(6, 2, 13);
			b7_sum += pc_2 * pow_mod(7, 2, 13);
			b11_sum += pc_2 * pow_mod(11, 2, 13);
			const size_t pc_3 = pop_count(upper_bits & (bitmask << 3));
			b6_sum += pc_3 * pow_mod(6, 3, 13);
			b7_sum += pc_3 * pow_mod(7, 3, 13);
			b11_sum += pc_3 * pow_mod(11, 3, 13);
			const size_t pc_4 = pop_count(upper_bits & (bitmask << 4));
			b6_sum += pc_4 * pow_mod(6, 4, 13);
			b7_sum += pc_4 * pow_mod(7, 4, 13);
			b11_sum += pc_4 * pow_mod(11, 4, 13);
			const size_t pc_5 = pop_count(upper_bits & (bitmask << 5));
			b6_sum += pc_5 * pow_mod(6, 5, 13);
			b7_sum += pc_5 * pow_mod(7, 5, 13);
			b11_sum += pc_5 * pow_mod(11, 5, 13);
			const size_t pc_6 = pop_count(upper_bits & (bitmask << 6));
			b6_sum += pc_6 * pow_mod(6, 6, 13);
			b7_sum += pc_6 * pow_mod(7, 6, 13);
			b11_sum += pc_6 * pow_mod(11, 6, 13);
			const size_t pc_7 = pop_count(upper_bits & (bitmask << 7));
			b6_sum += pc_7 * pow_mod(6, 7, 13);
			b7_sum += pc_7 * pow_mod(7, 7, 13);
			b11_sum += pc_7 * pow_mod(11, 7, 13);
			const size_t pc_8 = pop_count(upper_bits & (bitmask << 8));
			b6_sum += pc_8 * pow_mod(6, 8, 13);
			b7_sum += pc_8 * pow_mod(7, 8, 13);
			b11_sum += pc_8 * pow_mod(11, 8, 13);
			const size_t pc_9 = pop_count(upper_bits & (bitmask << 9));
			b6_sum += pc_9 * pow_mod(6, 9, 13);
			b7_sum += pc_9 * pow_mod(7, 9, 13);
			b11_sum += pc_9 * pow_mod(11, 9, 13);
			const size_t pc_10 = pop_count(upper_bits & (bitmask << 10));
			b6_sum += pc_10 * pow_mod(6, 10, 13);
			b7_sum += pc_10 * pow_mod(7, 10, 13);
			b11_sum += pc_10 * pow_mod(11, 10, 13);
			const size_t pc_11 = pop_count(upper_bits & (bitmask << 11));
			b6_sum += pc_11 * pow_mod(6, 11, 13);
			b7_sum += pc_11 * pow_mod(7, 11, 13);
			b11_sum += pc_11 * pow_mod(11, 11, 13);

			const uint8_t* const indivisible_b6 = indivisible_by[get_prime_index<13>::idx].data() + b6_sum;
			const uint8_t* const indivisible_b7 = indivisible_by[get_prime_index<13>::idx].data() + b7_sum;
			const uint8_t* const indivisible_b11 = indivisible_by[get_prime_index<13>::idx].data() + b11_sum;

			const size_t n_of_candidates = candidates_end - input;
			const size_t* const rounded_end = input + (n_of_candidates - (n_of_candidates % 2));

			const uint256_t shuffle_mask_lo = _mm256_set_epi64x(0x0909090909090909, 0x0808080808080808, 0x0101010101010101, 0x0000000000000000);
			const uint256_t shuffle_mask_hi = _mm256_set_epi64x(0x0B0B0B0B0B0B0B0B, 0x0A0A0A0A0A0A0A0A, 0x0303030303030303, 0x0202020202020202);
			const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

			const uint128_t xmm_b6_rems_lo = _mm_loadu_si128((uint128_t*)&static_b6_rems_lo);
			const uint256_t b6_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b6_rems_lo), xmm_b6_rems_lo, 1);
			const uint128_t xmm_b6_rems_hi = _mm_loadu_si128((uint128_t*)&static_b6_rems_hi);
			const uint256_t b6_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b6_rems_hi), xmm_b6_rems_hi, 1);
			const uint128_t xmm_b7_rems_lo = _mm_loadu_si128((uint128_t*)&static_b7_rems_lo);
			const uint256_t b7_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b7_rems_lo), xmm_b7_rems_lo, 1);
			const uint128_t xmm_b7_rems_hi = _mm_loadu_si128((uint128_t*)&static_b7_rems_hi);
			const uint256_t b7_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b7_rems_hi), xmm_b7_rems_hi, 1);
			const uint128_t xmm_b11_rems_lo = _mm_loadu_si128((uint128_t*)&static_b11_rems_lo);
			const uint256_t b11_rems_lo = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b11_rems_lo), xmm_b11_rems_lo, 1);
			const uint128_t xmm_b11_rems_hi = _mm_loadu_si128((uint128_t*)&static_b11_rems_hi);
			const uint256_t b11_rems_hi = _mm256_inserti128_si256(_mm256_castsi128_si256(xmm_b11_rems_hi), xmm_b11_rems_hi, 1);

			alignas(32) volatile uint32_t sums[8]{};

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

				// zero out garbage values so we can read results as u32s
				b6_7_11_x_sums = _mm256_blend_epi16(b6_7_11_x_sums, _mm256_setzero_si256(), 0b11101010);
				// sums are now stored as [0, 1, 2, x][0, 1, 2, x]

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

				b6_7_11_x_sums = _mm256_blend_epi16(b6_7_11_x_sums, _mm256_setzero_si256(), 0b11101010);

				// Only advance the pointer if the number is still a candidate

				const size_t inc_0 = indivisible_b6[sums[0]]
					& indivisible_b7[sums[1]]
					& indivisible_b11[sums[2]];
				*output = *input++;
				output = (uint64_t*)(((uint8_t*)output) + inc_0);

				const size_t inc_1 = indivisible_b6[sums[0 + 4]]
					& indivisible_b7[sums[1 + 4]]
					& indivisible_b11[sums[2 + 4]];
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
			b6_sum += pc_1 * pow_mod(6, 1, 13);
			b7_sum += pc_1 * pow_mod(7, 1, 13);
			b11_sum += pc_1 * pow_mod(11, 1, 13);
			const size_t pc_2 = pop_count(number & (bitmask << 2));
			b6_sum += pc_2 * pow_mod(6, 2, 13);
			b7_sum += pc_2 * pow_mod(7, 2, 13);
			b11_sum += pc_2 * pow_mod(11, 2, 13);
			const size_t pc_3 = pop_count(number & (bitmask << 3));
			b6_sum += pc_3 * pow_mod(6, 3, 13);
			b7_sum += pc_3 * pow_mod(7, 3, 13);
			b11_sum += pc_3 * pow_mod(11, 3, 13);
			const size_t pc_4 = pop_count(number & (bitmask << 4));
			b6_sum += pc_4 * pow_mod(6, 4, 13);
			b7_sum += pc_4 * pow_mod(7, 4, 13);
			b11_sum += pc_4 * pow_mod(11, 4, 13);
			const size_t pc_5 = pop_count(number & (bitmask << 5));
			b6_sum += pc_5 * pow_mod(6, 5, 13);
			b7_sum += pc_5 * pow_mod(7, 5, 13);
			b11_sum += pc_5 * pow_mod(11, 5, 13);
			const size_t pc_6 = pop_count(number & (bitmask << 6));
			b6_sum += pc_6 * pow_mod(6, 6, 13);
			b7_sum += pc_6 * pow_mod(7, 6, 13);
			b11_sum += pc_6 * pow_mod(11, 6, 13);
			const size_t pc_7 = pop_count(number & (bitmask << 7));
			b6_sum += pc_7 * pow_mod(6, 7, 13);
			b7_sum += pc_7 * pow_mod(7, 7, 13);
			b11_sum += pc_7 * pow_mod(11, 7, 13);
			const size_t pc_8 = pop_count(number & (bitmask << 8));
			b6_sum += pc_8 * pow_mod(6, 8, 13);
			b7_sum += pc_8 * pow_mod(7, 8, 13);
			b11_sum += pc_8 * pow_mod(11, 8, 13);
			const size_t pc_9 = pop_count(number & (bitmask << 9));
			b6_sum += pc_9 * pow_mod(6, 9, 13);
			b7_sum += pc_9 * pow_mod(7, 9, 13);
			b11_sum += pc_9 * pow_mod(11, 9, 13);
			const size_t pc_10 = pop_count(number & (bitmask << 10));
			b6_sum += pc_10 * pow_mod(6, 10, 13);
			b7_sum += pc_10 * pow_mod(7, 10, 13);
			b11_sum += pc_10 * pow_mod(11, 10, 13);
			const size_t pc_11 = pop_count(number & (bitmask << 11));
			b6_sum += pc_11 * pow_mod(6, 11, 13);
			b7_sum += pc_11 * pow_mod(7, 11, 13);
			b11_sum += pc_11 * pow_mod(11, 11, 13);

			// Only advance the pointer if the number is still a candidate
			const size_t inc = indivisible_by_13[b6_sum]
				& indivisible_by_13[b7_sum]
				& indivisible_by_13[b11_sum];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}

	template<bool on_fast_path>
	inline_toggle size_t* div_tests_with_16_rems(size_t* input,
													 const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t upper_bits_mask = size_t(-1) << 32;

		constexpr static std::array<uint8_t, 32> static_rems_01{
			1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6,     // base 3 % 17
			1, 5, 8, 6, 13, 14, 2, 10, 16, 12, 9, 11, 4, 3, 15, 7 }; // base 5 % 17
		constexpr static std::array<uint8_t, 32> static_rems_23{
			1, 6, 2, 12, 4, 7, 8, 14, 16, 11, 15, 5, 13, 10, 9, 3,     // base 6 % 17
			1, 7, 15, 3, 4, 11, 9, 12, 16, 10, 2, 14, 13, 6, 8, 5 }; // base 7 % 17
		constexpr static std::array<uint8_t, 32> static_rems_45{
			1, 10, 15, 14, 4, 6, 9, 5, 16, 7, 2, 3, 13, 11, 8, 12,     // base 10 % 17
			1, 11, 2, 5, 4, 10, 8, 3, 16, 6, 15, 12, 13, 7, 9, 14 }; // base 11 % 17
		constexpr static std::array<uint8_t, 32> static_rems_6x{
			1, 12, 8, 11, 13, 3, 2, 7, 16, 5, 9, 6, 4, 14, 15, 10,     // base 12 % 17
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // padding

		const uint8_t* const indivisible_by_17 = indivisible_by[get_prime_index<17>::idx].data();

		size_t* output = input;
		size_t number = *input;

		size_t upper_bits = number & upper_bits_mask;
		// get popcounts of upper 32 bits, in both lanes, in range 0-2
		uint256_t upper_bytes = util::expand_bits_to_bytes(number >> 32);
		upper_bytes = _mm256_and_si256(upper_bytes, _mm256_set1_epi8(0x01)); // convert 0xFF -> 0x01
		upper_bytes = _mm256_add_epi8(upper_bytes, _mm256_permute2x128_si256(upper_bytes, upper_bytes, 1));

		const uint256_t rems_01 = _mm256_loadu_si256((uint256_t*)&static_rems_01);
		const uint256_t rems_23 = _mm256_loadu_si256((uint256_t*)&static_rems_23);
		const uint256_t rems_45 = _mm256_loadu_si256((uint256_t*)&static_rems_45);
		const uint256_t rems_6x = _mm256_loadu_si256((uint256_t*)&static_rems_6x);

		const uint256_t shuffle_mask = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, 0x0101010101010101, 0x0000000000000000);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		alignas(32) volatile uint32_t sums[8]{};

		{
			// load one candidate into all positions
			const uint256_t candidate = _mm256_set1_epi64x(*input);

			// convert bits to bytes
			uint256_t lower_bytes = _mm256_shuffle_epi8(candidate, shuffle_mask);
			lower_bytes = _mm256_andnot_si256(lower_bytes, and_mask);
			lower_bytes = _mm256_cmpeq_epi8(lower_bytes, _mm256_setzero_si256());

			// get popcounts of all 64 bits in range 0-4, in both lanes
			uint256_t pc = _mm256_sub_epi8(upper_bytes, lower_bytes);
			pc = _mm256_sub_epi8(pc, _mm256_permute2x128_si256(lower_bytes, lower_bytes, 1));

			// multiply 16+16 8-bit values by 16+16 remainders, storing partially summed results as 8+8 uint16_ts
			const uint256_t sums_01 = _mm256_maddubs_epi16(pc, rems_01);
			const uint256_t sums_23 = _mm256_maddubs_epi16(pc, rems_23);
			const uint256_t sums_45 = _mm256_maddubs_epi16(pc, rems_45);
			const uint256_t sums_6x = _mm256_maddubs_epi16(pc, rems_6x);

			// pack into 4x 8-bit integers
			uint256_t sums_0213 = _mm256_packus_epi16(sums_01, sums_23);
			uint256_t sums_465x = _mm256_packus_epi16(sums_45, sums_6x);

			// h-sum into 4x 16-bit integers, 2 in each lane
			sums_0213 = _mm256_sad_epu8(sums_0213, _mm256_setzero_si256());
			sums_465x = _mm256_sad_epu8(sums_465x, _mm256_setzero_si256());

			// pack to 8x 32-bit integers
			const uint256_t sums_0246135x = _mm256_packus_epi32(sums_0213, sums_465x);

			// store on the stack
			_mm256_storeu_si256((uint256_t*)sums, sums_0246135x);
		}

		// load ahead
		uint256_t candidate = _mm256_set1_epi64x(*(input + 1));

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

			uint256_t lower_bytes = _mm256_shuffle_epi8(candidate, shuffle_mask);

			// load two iterations ahead
			candidate = _mm256_set1_epi64x(*(input + 2));

			lower_bytes = _mm256_andnot_si256(lower_bytes, and_mask);
			lower_bytes = _mm256_cmpeq_epi8(lower_bytes, _mm256_setzero_si256());

			uint256_t pc = _mm256_sub_epi8(upper_bytes, lower_bytes);
			pc = _mm256_sub_epi8(pc, _mm256_permute2x128_si256(lower_bytes, lower_bytes, 1));

			const uint256_t sums_01 = _mm256_maddubs_epi16(pc, rems_01);
			const uint256_t sums_23 = _mm256_maddubs_epi16(pc, rems_23);
			const uint256_t sums_45 = _mm256_maddubs_epi16(pc, rems_45);
			const uint256_t sums_6x = _mm256_maddubs_epi16(pc, rems_6x);

			uint256_t sums_0213 = _mm256_packus_epi16(sums_01, sums_23);
			uint256_t sums_465x = _mm256_packus_epi16(sums_45, sums_6x);

			sums_0213 = _mm256_sad_epu8(sums_0213, _mm256_setzero_si256());
			sums_465x = _mm256_sad_epu8(sums_465x, _mm256_setzero_si256());

			const uint256_t sums_0246135x = _mm256_packus_epi32(sums_0213, sums_465x);

			*output = *input++; // always write

			// Only advance the pointer if the number is still a candidate
			const size_t inc = indivisible_by_17[sums[0]]
				& indivisible_by_17[sums[1]]
				& indivisible_by_17[sums[2]]
				& indivisible_by_17[sums[3]]
				& indivisible_by_17[sums[4]]
				& indivisible_by_17[sums[5]]
				& indivisible_by_17[sums[6]];
			output = (uint64_t*)(((uint8_t*)output) + inc);

			// store the above results for the next iteration
			_mm256_storeu_si256((uint256_t*)sums, sums_0246135x);
		}

		return output;
	}

	template<bool on_fast_path>
	inline_toggle size_t* div_tests_with_8_rems(size_t* input,
													const size_t* const candidates_end)
	{
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t bitmask = bitmask_for<8, 17>::val;
		static_assert(bitmask == bitmask_for<9, 17>::val);
		static_assert(period_of<bitmask>() == 8);

		size_t* output = input;

		if constexpr (on_fast_path)
		{
			constexpr static std::array<uint8_t, 32> static_nybble_lookup_8m17_lo = mbp::detail::build_8rem_shuffle_lookup_lo_nybble<8, 17>();
			constexpr static std::array<uint8_t, 32> static_nybble_lookup_8m17_hi = mbp::detail::build_8rem_shuffle_lookup_hi_nybble<8, 17>();
			constexpr static std::array<uint8_t, 32> static_nybble_lookup_9m17_lo = mbp::detail::build_8rem_shuffle_lookup_lo_nybble<9, 17>();
			constexpr static std::array<uint8_t, 32> static_nybble_lookup_9m17_hi = mbp::detail::build_8rem_shuffle_lookup_hi_nybble<9, 17>();

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_lookup_8m17_lo = _mm256_loadu_si256((uint256_t*)&static_nybble_lookup_8m17_lo);
			const uint256_t nybble_lookup_8m17_hi = _mm256_loadu_si256((uint256_t*)&static_nybble_lookup_8m17_hi);
			const uint256_t nybble_lookup_9m17_lo = _mm256_loadu_si256((uint256_t*)&static_nybble_lookup_9m17_lo);
			const uint256_t nybble_lookup_9m17_hi = _mm256_loadu_si256((uint256_t*)&static_nybble_lookup_9m17_hi);

			const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) % 8);

			// the upper 32 bits don't change, calculate their sum of remainders
			const size_t number_upper = (*input) & (size_t(-1) << 32);
			const size_t pc_0 = pop_count(number_upper & (bitmask << 0));
			size_t b8m17_rem = pc_0;
			size_t b9m17_rem = pc_0;
			const size_t pc_1 = pop_count(number_upper & (bitmask << 1));
			b8m17_rem += pc_1 * pow_mod(8, 1, 17);
			b9m17_rem += pc_1 * pow_mod(9, 1, 17);
			const size_t pc_2 = pop_count(number_upper & (bitmask << 2));
			b8m17_rem += pc_2 * pow_mod(8, 2, 17);
			b9m17_rem += pc_2 * pow_mod(9, 2, 17);
			const size_t pc_3 = pop_count(number_upper & (bitmask << 3));
			b8m17_rem += pc_3 * pow_mod(8, 3, 17);
			b9m17_rem += pc_3 * pow_mod(9, 3, 17);
			const size_t pc_4 = pop_count(number_upper & (bitmask << 4));
			b8m17_rem += pc_4 * pow_mod(8, 4, 17);
			b9m17_rem += pc_4 * pow_mod(9, 4, 17);
			const size_t pc_5 = pop_count(number_upper & (bitmask << 5));
			b8m17_rem += pc_5 * pow_mod(8, 5, 17);
			b9m17_rem += pc_5 * pow_mod(9, 5, 17);
			const size_t pc_6 = pop_count(number_upper & (bitmask << 6));
			b8m17_rem += pc_6 * pow_mod(8, 6, 17);
			b9m17_rem += pc_6 * pow_mod(9, 6, 17);
			const size_t pc_7 = pop_count(number_upper & (bitmask << 7));
			b8m17_rem += pc_7 * pow_mod(8, 7, 17);
			b9m17_rem += pc_7 * pow_mod(9, 7, 17);
			const uint8_t* const indivisible_b8 = indivisible_by[get_prime_index<17>::idx].data() + b8m17_rem;
			const uint8_t* const indivisible_b9 = indivisible_by[get_prime_index<17>::idx].data() + b9m17_rem;

			alignas(32) volatile uint32_t sums[16]{};

			// load 8 candidates into two registers
			uint256_t ymm0 = _mm256_loadu_si256((uint256_t*)(input + 0));
			uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 4));
			// keep the bottom 32 bits of our 8 candidates
			ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

			// run vector instructions one iteration ahead
			{
				// select the high and low halves of each byte
				uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// load two iterations ahead
				ymm0 = _mm256_loadu_si256((uint256_t*)(input + 8));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 12));
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

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

				// pack u64 results -> u32s
				const uint256_t sums_04152637_b8m17 = _mm256_packus_epi32(sums_0426_b8m17, sums_1537_b8m17);
				const uint256_t sums_04152637_b9m17 = _mm256_packus_epi32(sums_0426_b9m17, sums_1537_b9m17);

				// store results on the stack
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_04152637_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_04152637_b9m17);
			}

			for (; input < candidates_end_rounded; )
			{
				// run vector instructions one iteration ahead
				uint256_t lo_nybbles = _mm256_and_si256(ymm0, nybble_mask);
				uint256_t hi_nybbles = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

				// load two iterations ahead
				ymm0 = _mm256_loadu_si256((uint256_t*)(input + 16));
				uint256_t ymm1 = _mm256_loadu_si256((uint256_t*)(input + 20));
				ymm0 = _mm256_blend_epi32(ymm0, _mm256_slli_epi64(ymm1, 32), 0b10'10'10'10);

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

				const uint256_t sums_04152637_b8m17 = _mm256_packus_epi32(sums_0426_b8m17, sums_1537_b8m17);
				const uint256_t sums_04152637_b9m17 = _mm256_packus_epi32(sums_0426_b9m17, sums_1537_b9m17);

				*output = *input++; // always copy
				const size_t inc_a = indivisible_b8[sums[0]] & indivisible_b9[sums[0 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_a); // branchless conditional increment

				*output = *input++;
				const size_t inc_b = indivisible_b8[sums[2]] & indivisible_b9[sums[2 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_b);

				*output = *input++;
				const size_t inc_c = indivisible_b8[sums[4]] & indivisible_b9[sums[4 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_c);

				*output = *input++;
				const size_t inc_d = indivisible_b8[sums[6]] & indivisible_b9[sums[6 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_d);

				*output = *input++;
				const size_t inc_e = indivisible_b8[sums[1]] & indivisible_b9[sums[1 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_e);

				*output = *input++;
				const size_t inc_f = indivisible_b8[sums[3]] & indivisible_b9[sums[3 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_f);

				*output = *input++;
				const size_t inc_g = indivisible_b8[sums[5]] & indivisible_b9[sums[5 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_g);

				*output = *input++;
				const size_t inc_h = indivisible_b8[sums[7]] & indivisible_b9[sums[7 + 8]];
				output = (uint64_t*)(((uint8_t*)output) + inc_h);

				// store the above results for the next iteration
				_mm256_storeu_si256((uint256_t*)(sums + 0), sums_04152637_b8m17);
				_mm256_storeu_si256((uint256_t*)(sums + 8), sums_04152637_b9m17);
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
			b8m17_rem += pc_1 * pow_mod(8, 1, 17);
			b9m17_rem += pc_1 * pow_mod(9, 1, 17);
			const size_t pc_2 = pop_count(candidate & (bitmask << 2));
			b8m17_rem += pc_2 * pow_mod(8, 2, 17);
			b9m17_rem += pc_2 * pow_mod(9, 2, 17);
			const size_t pc_3 = pop_count(candidate & (bitmask << 3));
			b8m17_rem += pc_3 * pow_mod(8, 3, 17);
			b9m17_rem += pc_3 * pow_mod(9, 3, 17);
			const size_t pc_4 = pop_count(candidate & (bitmask << 4));
			b8m17_rem += pc_4 * pow_mod(8, 4, 17);
			b9m17_rem += pc_4 * pow_mod(9, 4, 17);
			const size_t pc_5 = pop_count(candidate & (bitmask << 5));
			b8m17_rem += pc_5 * pow_mod(8, 5, 17);
			b9m17_rem += pc_5 * pow_mod(9, 5, 17);
			const size_t pc_6 = pop_count(candidate & (bitmask << 6));
			b8m17_rem += pc_6 * pow_mod(8, 6, 17);
			b9m17_rem += pc_6 * pow_mod(9, 6, 17);
			const size_t pc_7 = pop_count(candidate & (bitmask << 7));
			b8m17_rem += pc_7 * pow_mod(8, 7, 17);
			b9m17_rem += pc_7 * pow_mod(9, 7, 17);

			// conditionally increment
			const auto inc = indivisible_by_17[b8m17_rem] & indivisible_by_17[b9m17_rem];
			output = (uint64_t*)(((uint8_t*)output) + inc);
		}

		return output;
	}
}
