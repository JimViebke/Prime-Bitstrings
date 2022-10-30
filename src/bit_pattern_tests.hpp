#pragma once

#include "config.hpp"
#include "math/math.hpp"

namespace mbp
{
	__forceinline void merge_static_sieve_with_bit_pattern_filters(
		uint256_t* const out,
		const uint256_t* const static_sieve_in,
		const uint256_t* const pc_in,
		const uint256_t& valid_pcs,
		const uint256_t* const gcd_in)
	{
		uint256_t pc_mask_0 = _mm256_load_si256(pc_in + 0);
		uint256_t ymm0 = _mm256_load_si256(static_sieve_in + 0);
		uint256_t gcd_mask_0 = _mm256_loadu_si256(gcd_in + 0);
		pc_mask_0 = _mm256_shuffle_epi8(valid_pcs, pc_mask_0);
		ymm0 = _mm256_and_si256(ymm0, pc_mask_0);
		ymm0 = _mm256_and_si256(ymm0, gcd_mask_0);
		_mm256_store_si256(out + 0, ymm0);

		uint256_t pc_mask_1 = _mm256_load_si256(pc_in + 1);
		uint256_t ymm1 = _mm256_load_si256(static_sieve_in + 1);
		uint256_t gcd_mask_1 = _mm256_loadu_si256(gcd_in + 1);
		pc_mask_1 = _mm256_shuffle_epi8(valid_pcs, pc_mask_1);
		ymm1 = _mm256_and_si256(ymm1, pc_mask_1);
		ymm1 = _mm256_and_si256(ymm1, gcd_mask_1);
		_mm256_store_si256(out + 1, ymm1);

		uint256_t pc_mask_2 = _mm256_load_si256(pc_in + 2);
		uint256_t ymm2 = _mm256_load_si256(static_sieve_in + 2);
		uint256_t gcd_mask_2 = _mm256_loadu_si256(gcd_in + 2);
		pc_mask_2 = _mm256_shuffle_epi8(valid_pcs, pc_mask_2);
		ymm2 = _mm256_and_si256(ymm2, pc_mask_2);
		ymm2 = _mm256_and_si256(ymm2, gcd_mask_2);
		_mm256_store_si256(out + 2, ymm2);

		uint256_t pc_mask_3 = _mm256_load_si256(pc_in + 3);
		uint256_t ymm3 = _mm256_load_si256(static_sieve_in + 3);
		uint256_t gcd_mask_3 = _mm256_loadu_si256(gcd_in + 3);
		pc_mask_3 = _mm256_shuffle_epi8(valid_pcs, pc_mask_3);
		ymm3 = _mm256_and_si256(ymm3, pc_mask_3);
		ymm3 = _mm256_and_si256(ymm3, gcd_mask_3);
		_mm256_store_si256(out + 3, ymm3);
	}

	uint256_t manually_generate_bit_pattern_filters(size_t number)
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();
		constexpr size_t tiny_gcd_lookup = []() consteval {
			size_t val = 0;
			for (size_t i = 0; i < 32; ++i)
				val |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << i;
			return val;
		}();

		uint256_t mask{};

		for (size_t idx = 0; idx < 32; ++idx, number += 2)
		{
			const auto valid_pc = tiny_primes_lookup >> pop_count(number); // & 1

			const auto valid_gcd = tiny_gcd_lookup >> abs(pop_count(number & 0x5555555555555555) -
														  pop_count(number & 0xAAAAAAAAAAAAAAAA)); // & 1

			mask.m256i_u8[idx] = valid_pc & valid_gcd & 1u;
		}

		return mask;
	}



	std::vector<uint8_t> build_popcounts_lookup()
	{
		std::vector<uint8_t> popcounts;
		popcounts.reserve(1ull << 16);

		for (size_t i = 0; i < (1ull << 16); ++i)
		{
			popcounts.push_back(uint8_t(pop_count(i)));
		}

		return popcounts;
	}
	const std::vector<uint8_t> pc_lookup = build_popcounts_lookup();

	std::vector<std::vector<uint8_t>> build_gcd_lookup()
	{
		/*
		We have 48 "outer" bits (47 upper and 1 lower) and 16 "inner" bits.

		With 48 outer bits, we have 24 even and 24 odd bits:
			even bits can be 1 through 24	(bottom bit is always set)
			odd  bits can be 0 through 24
			Added together, (even - odd) yields range: -23 to +24

		Normalize outer alternating bitsum from -23,+24 to 0,47 by adding 23

		Build 48 lookup tables, each containing 2^16 bytes. The byte at index i reads 1 if the bit pattern of
		i allows for a valid GCD, 0 otherwise.
		*/

		// contains 1 bits at indexes that share a GCD of 1 with a product of primes
		constexpr size_t tiny_gcd_lookup = []() consteval {
			size_t val = 0;
			for (size_t i = 0; i < 32; ++i)
			{
				val |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << i;
			}
			return val;
		}();

		std::vector<std::vector<uint8_t>> gcd_lookup;
		gcd_lookup.reserve(48);

		for (int outer_abs = -23; outer_abs <= 24; ++outer_abs)
		{
			std::vector<uint8_t> lookup;
			lookup.reserve(1ull << 16);

			for (size_t i = 0; i < (1ull << 16); ++i)
			{
				// calculate the alternating bitsum of i, add "outer_abs"
				const auto even_pc = pop_count(i & 0xAAAAAAAAAAAAAAAA); // We're generating a lookup for bits 1-16,
				const auto odd_pc = pop_count(i & 0x5555555555555555); //  so the even/odd masks are swapped
				const auto alternating_bitsum = (even_pc - odd_pc) + outer_abs;

				lookup.push_back(uint8_t((tiny_gcd_lookup >> abs(alternating_bitsum)) & 1));
			}

			gcd_lookup.push_back(lookup);
		}

		return gcd_lookup;
	}
	const std::vector<std::vector<uint8_t>> gcd_lookup = build_gcd_lookup();

}
