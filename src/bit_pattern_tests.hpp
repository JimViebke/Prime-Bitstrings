#pragma once

#include "config.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"

namespace mbp
{
	constexpr size_t pow_2_16 = 1ull << 16;

	std::vector<std::vector<bit_array<pow_2_16>>> build_popcounts_lookup()
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();

		std::vector<std::vector<bit_array<pow_2_16>>> pc_lookup;
		pc_lookup.reserve(8);

		// for each offset
		for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset)
		{
			auto& offset_lookup = pc_lookup.emplace_back();
			// subtract 2 to normalize outer popcount of 2-48 to idx 0-46 (47 elements)
			offset_lookup.reserve(47);

			// for each outer popcount
			for (size_t outer_pc = 2; outer_pc <= 48; ++outer_pc)
			{
				const size_t shifted_primes_lookup = tiny_primes_lookup >> outer_pc;
				uint64_t* bit_array_ptr = (uint64_t*)offset_lookup.emplace_back().data();

				// inner bits == bits 1 through 16 (not 0 through 15)

				// Start inner_bits at 0-7 instead of 0. This drops off the first 0-7 (unused) values,
				// and has the same effect as shifting this part of the lookup left by 0-7 bits.
				for (size_t inner_bits = bit_offset; inner_bits < pow_2_16; /* increment below */)
				{
					uint64_t chunk = 0;
					for (size_t j = 0; j < 64; ++j, ++inner_bits)
					{
						const size_t bit = (shifted_primes_lookup >> pop_count(inner_bits)) & 1;
						chunk |= (bit << j);
					}
					*bit_array_ptr = chunk;
					++bit_array_ptr;
				}
			}
		}

		return pc_lookup;
	}
	// [offset 0-7][outer pc 0-46][bitstring inner bits]
	const std::vector<std::vector<bit_array<pow_2_16>>> pc_lookup = build_popcounts_lookup();

	std::vector<std::vector<bit_array<pow_2_16>>> build_gcd_lookup()
	{
		/*
		We have 48 "outer" bits (47 upper and 1 lower) and 16 "inner" bits (bits 1 through 16; not 0 through 15)

		With 48 outer bits, we have 24 even and 24 odd bits:
			even bits can be 1 through 24	(bottom bit is always set)
			odd  bits can be 0 through 24
			Added together, (even - odd) yields range: -23 to +24

		Normalize outer alternating bitsum from -23,+24 to 0,47 by adding 23
		*/

		// contains a 1 bit at indexes that share a GCD of 1 with a product of primes
		constexpr size_t tiny_gcd_lookup = []() consteval {
			size_t val = 0;
			for (size_t i = 0; i < 32; ++i)
			{
				val |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << i;
			}
			return val;
		}();

		std::vector<std::vector<bit_array<pow_2_16>>> gcd_lookup;
		gcd_lookup.reserve(8);

		// for each offset
		for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset)
		{
			auto& offset_lookup = gcd_lookup.emplace_back();
			offset_lookup.reserve(48);

			// for each outer alternating bitsum
			for (int outer_abs = -23; outer_abs <= 24; ++outer_abs)
			{
				uint64_t* bit_array_ptr = (uint64_t*)offset_lookup.emplace_back().data();

				// Start inner_bits at 0-7 instead of 0. This drops off the first 0-7 (unused) values,
				// and has the same effect as shifting this part of the lookup left by 0-7 bits.
				for (size_t inner_bits = bit_offset; inner_bits < pow_2_16; /* increment below */)
				{
					// combine the next 64 writes
					uint64_t chunk = 0;
					for (size_t j = 0; j < 64; ++j, ++inner_bits)
					{
						// calculate the alternating bitsum of the inner bits, add outer_abs
						const auto even_pc = pop_count(inner_bits & 0xAAAAAAAAAAAAAAAA); // We're generating a lookup for bits 1-16,
						const auto odd_pc = pop_count(inner_bits & 0x5555555555555555); //  so the even/odd masks are swapped
						const auto alternating_bitsum = (even_pc - odd_pc) + outer_abs;

						const size_t bit = (tiny_gcd_lookup >> abs(alternating_bitsum)) & 1;
						chunk |= (bit << j);
					}

					*bit_array_ptr = chunk;
					++bit_array_ptr;
				}
			}
		}

		return gcd_lookup;
	}
	// [offset 0-7][outer abs 0-47][bitstring inner bits]
	const std::vector<std::vector<bit_array<pow_2_16>>> gcd_lookup = build_gcd_lookup();

	std::array<std::vector<bit_array<pow_2_16>>, 8> build_b3m5_lookup()
	{
		std::array<std::vector<bit_array<pow_2_16>>, 8> lookup{};

		// 1, 3, 4, 2, 1, 3, 4, 2, 1, 3, 4, 2, 1, 3, 4, 2, 1
		//    ^ start here
		constexpr std::array<uint8_t, 16> nybble_lookup = []() consteval {
			std::array<uint8_t, 16> arr{};
			for (size_t i = 0; i < 16; ++i)
			{
				if (i & 0b0001) arr[i] += 3;
				if (i & 0b0010) arr[i] += 4;
				if (i & 0b0100) arr[i] += 2;
				if (i & 0b1000) arr[i] += 1;
			}
			return arr;
		}();

		for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset)
		{
			auto& offset_lookup = lookup[bit_offset];
			offset_lookup.reserve(5);

			for (size_t outer_rem = 0; outer_rem < 5; ++outer_rem)
			{
				auto& bit_array = offset_lookup.emplace_back();

				// Start at 0-7 instead of 0. This drops off the first 0-7 (unused) values,
				// and has the same effect as shifting this part of the lookup left by 0-7 bits.
				// This also leave the upper 0-7 bits unset, which we will never read because they
				// are past the 2^16 rollover.
				auto nybble_d_idx = bit_offset;
				size_t i = 0;

				// 16^4 == the 65,536 bits we're looking for
				for (auto nybble_a : nybble_lookup)
				{
					const auto sum_a = outer_rem + nybble_a;

					for (auto nybble_b : nybble_lookup)
					{
						const auto sum_ab = sum_a + nybble_b;

						for (auto nybble_c : nybble_lookup)
						{
							const auto sum_abc = sum_ab + nybble_c;

							// nybble_d_idx starts at [bit_offset], then 0 for all other iterations
							for (; nybble_d_idx < 16; ++nybble_d_idx)
							{
								// largest possible sum of four nybbles + outer_rem == 44 (45-bit lookup table)
								constexpr size_t indivisible_by_5 = 0b11110'11110'11110'11110'11110'11110'11110'11110'11110;

								// calculate the remainder sum of all 16 inner bits + the outer bits
								const auto sum_abcd = sum_abc + nybble_lookup[nybble_d_idx];
								// extract a bit from our lookup
								const auto indivisible = (indivisible_by_5 >> sum_abcd) & 1;

								// If a candidate is indivisible by 5 in base 3, it is still a candidate
								// Set that bit high.
								bit_array.data()[i / 8] |= (indivisible << (i % 8));

								++i;
							}
							nybble_d_idx = 0;
						}
					}
				}
			}
		}

		return lookup;
	}
	// [offset 0-7][outer rem 0-4][bitstring inner bits]
	const std::array<std::vector<bit_array<pow_2_16>>, 8> b3m5_lookup = build_b3m5_lookup();

}
